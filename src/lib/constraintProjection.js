// src/lib/constraintProjection.js
/**
 * Implements constraint projection using the saddle point system
 * Following Sections 5.3.1 and 5.3.2 of the paper
 */

import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from './stores';
import { build_A_bar_2D } from './innerProduct';

// Numerical stability parameters
export const ProjectionNumerics = {
    // Toggle numerical stability measures
    USE_REGULARIZATION: false,         // Use regularization in saddle point system
    USE_ROBUST_SOLVER: false,          // Use more robust solver techniques
    LOG_MATRIX_PROPERTIES: true,      // Log matrix properties for debugging
    
    // Parameter values
    SADDLE_REGULARIZATION: 1e-8,      // Regularization term for saddle point system
    MIN_PIVOT: 1e-12,                 // Minimum pivot value for direct solver
    MAX_CORRECTION_NORM: 1000.0,      // Maximum correction norm to prevent explosions
    DAMPING_FACTOR: 0.5,              // Damping factor for constraint projection
    CONSTRAINT_TOL: 1e-6              // Tolerance for constraint satisfaction
};

/**
 * Check the conditioning of a matrix
 * Helps diagnose numerical issues
 * 
 * @param {Array} matrix - Matrix to check
 * @param {string} name - Matrix name for logging
 */
function checkMatrixConditioning(matrix, name) {
    if (!ProjectionNumerics.LOG_MATRIX_PROPERTIES) return;
    
    console.log(`==== CHECKING MATRIX CONDITIONING: ${name} ====`);
    
    try {
        const mathMatrix = math.matrix(matrix);
        const size = mathMatrix.size();
        const rows = size[0];
        const cols = size.length > 1 ? size[1] : 1;
        
        console.log(`Matrix dimensions: ${rows} × ${cols}`);
        
        // Find min, max, and mean of absolute values
        let min = Infinity;
        let max = 0;
        let sum = 0;
        let count = 0;
        
        mathMatrix.forEach((value) => {
            const absValue = Math.abs(value);
            if (absValue > 0) {
                min = Math.min(min, absValue);
                max = Math.max(max, absValue);
                sum += absValue;
                count++;
            }
        });
        
        const mean = count > 0 ? sum / count : 0;
        console.log(`Non-zero values: min=${min.toExponential(4)}, max=${max.toExponential(4)}, mean=${mean.toExponential(4)}, count=${count}`);
        
        // Try to estimate condition number for square matrices
        if (rows === cols && rows <= 100) { // Limit to small matrices for performance
            try {
                // This is expensive and might fail for large or ill-conditioned matrices
                const eigs = math.eigs(mathMatrix);
                const eigenvalues = eigs.values;
                
                // Get real parts and filter non-zero values
                const realEigs = eigenvalues.filter(v => Math.abs(v) > 1e-14);
                
                if (realEigs.length > 0) {
                    const minEig = Math.min(...realEigs.map(Math.abs));
                    const maxEig = Math.max(...realEigs.map(Math.abs));
                    const condEst = maxEig / minEig;
                    console.log(`Eigenvalue-based condition number estimate: ${condEst.toExponential(4)}`);
                    
                    if (condEst > 1e10) {
                        console.warn(`Matrix ${name} is likely ill-conditioned!`);
                    }
                }
            } catch (e) {
                console.log(`Could not compute eigenvalues: ${e.message}`);
            }
        } else if (rows === cols) {
            // For larger matrices, check diagonal dominance as a rough heuristic
            let isDiagDominant = true;
            
            for (let i = 0; i < rows; i++) {
                const diagValue = Math.abs(matrix[i][i]);
                let rowSum = 0;
                
                for (let j = 0; j < cols; j++) {
                    if (i !== j) rowSum += Math.abs(matrix[i][j]);
                }
                
                if (diagValue <= rowSum) {
                    isDiagDominant = false;
                    break;
                }
            }
            
            console.log(`Matrix is ${isDiagDominant ? '' : 'not '}diagonally dominant`);
        }
        
        // Check for symmetry in square matrices
        if (rows === cols) {
            let isSymmetric = true;
            let maxAsymmetry = 0;
            
            for (let i = 0; i < rows && isSymmetric; i++) {
                for (let j = i+1; j < cols; j++) {
                    const diff = Math.abs(matrix[i][j] - matrix[j][i]);
                    maxAsymmetry = Math.max(maxAsymmetry, diff);
                    if (diff > 1e-10) {
                        isSymmetric = false;
                        break;
                    }
                }
            }
            
            console.log(`Matrix is ${isSymmetric ? '' : 'not '}symmetric, max asymmetry: ${maxAsymmetry.toExponential(4)}`);
        }
    } catch (error) {
        console.error(`Error checking matrix conditioning: ${error.message}`);
    }
    
    console.log("======================================");
}

/**
 * Solves the saddle point system for constrained optimization
 * [ A C^T ] [ x ] = [ b ]
 * [ C  0  ] [ λ ] = [ d ]
 *
 * @param {Array} A - Matrix A (inner product matrix)
 * @param {Array} C - Constraint Jacobian
 * @param {Array} b - RHS vector b 
 * @param {Array} d - RHS vector d (constraint values)
 * @returns {Object} - Solution vectors x and λ
 */
export function solveSaddlePointSystem(A, C, b, d) {
    console.log("==== SOLVING SADDLE POINT SYSTEM ====");
    console.log(`A dimensions: ${A.length} × ${A[0]?.length || 0}`);
    console.log(`C dimensions: ${C.length} × ${C[0]?.length || 0}`);
    console.log(`b length: ${b.length}, d length: ${d.length}`);
    
    try {
        const n = A.length;          // Dimension of primal variables
        const m = C.length;          // Number of constraints
        
        if (m === 0) {
            console.log("No constraints, returning solution of Ax = b");
            const x = math.lusolve(A, b);
            console.log("================================");
            return {
                x: math.flatten(x)._data,
                lambda: []
            };
        }
        
        console.log("Setting up saddle point system");
        
        /*
         * The saddle point system has the form:
         * [ A   C^T ] [ x ] = [ b ]
         * [ C    0  ] [ λ ] = [ d ]
         * 
         * For improved numerical stability, we can add a small regularization term:
         * [ A   C^T ] [ x ] = [ b ]
         * [ C   -εI ] [ λ ] = [ d ]
         */
        
        // Add regularization if enabled
        const epsilon = ProjectionNumerics.USE_REGULARIZATION ? 
            ProjectionNumerics.SADDLE_REGULARIZATION : 0;
        
        if (epsilon > 0) {
            console.log(`Using regularization with ε = ${epsilon}`);
        }
        
        // Build the full saddle point matrix
        const fullSize = n + m;
        const fullMatrix = [];
        
        // Initialize full matrix with zeros
        for (let i = 0; i < fullSize; i++) {
            fullMatrix[i] = new Array(fullSize).fill(0);
        }
        
        // Fill A block
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                fullMatrix[i][j] = A[i][j];
            }
        }
        
        // Fill C^T block
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < m; j++) {
                fullMatrix[i][n + j] = C[j][i]; // Transpose
            }
        }
        
        // Fill C block
        for (let i = 0; i < m; i++) {
            for (let j = 0; j < n; j++) {
                fullMatrix[n + i][j] = C[i][j];
            }
        }
        
        // Fill zeros block with regularization
        for (let i = 0; i < m; i++) {
            fullMatrix[n + i][n + i] = -epsilon;
        }
        
        // Build full right-hand side vector
        const fullRHS = b.concat(d);
        
        // Check matrix conditioning
        checkMatrixConditioning(fullMatrix, "Saddle Point Matrix");
        
        // Solve the system
        console.log("Solving saddle point system using math.js...");
        
        let solution;
        try {
            // Try direct solver with LU decomposition
            solution = math.lusolve(fullMatrix, fullRHS);
            console.log("Solved using LU decomposition");
        } catch (error) {
            if (ProjectionNumerics.USE_ROBUST_SOLVER) {
                console.warn(`LU solver failed: ${error.message}, trying fallback solver`);
                
                // Add small value to diagonal for stability if needed
                for (let i = 0; i < fullSize; i++) {
                    if (Math.abs(fullMatrix[i][i]) < ProjectionNumerics.MIN_PIVOT) {
                        const oldValue = fullMatrix[i][i];
                        fullMatrix[i][i] = Math.sign(oldValue || 1) * ProjectionNumerics.MIN_PIVOT;
                        console.log(`Reinforced diagonal at (${i},${i}): ${oldValue} -> ${fullMatrix[i][i]}`);
                    }
                }
                
                // Try again with modified system
                solution = math.lusolve(fullMatrix, fullRHS);
                console.log("Solved using LU decomposition with reinforced diagonal");
            } else {
                throw error;
            }
        }
        
        let x, lambda;
        try {
            // Solution from math.lusolve is a Matrix, convert to array properly
            const solutionArray = math.matrix(solution).toArray();
            
            // Extract x and lambda components
            x = solutionArray.slice(0, n).map(row => row[0]); // Get first column of each row
            lambda = solutionArray.slice(n, n + m).map(row => row[0]);
            
            // Calculate norms safely
            const xNorm = math.norm(x);
            const lambdaNorm = math.norm(lambda);
            
            console.log(`Solution found: |x| = ${xNorm.toExponential(4)}, |λ| = ${lambdaNorm.toExponential(4)}`);
            
            // Check for unreasonably large values
            if (xNorm > 1e10 || lambdaNorm > 1e10) {
                console.warn("Solution has extremely large values, possible numerical issues");
            }
        } catch (normError) {
            console.warn(`Could not calculate solution norms: ${normError.message}`);
            console.log("Solution extraction failed, trying fallback method");
            
            // Fallback extraction method
            x = [];
            for (let i = 0; i < n; i++) {
                x.push(solution[i][0] || 0);
            }
            
            lambda = [];
            for (let i = 0; i < m; i++) {
                lambda.push(solution[n+i][0] || 0);
            }
        }
        
        console.log("======================================");
        return { x, lambda };
    } catch (error) {
        console.error(`Error solving saddle point system: ${error.message}`);
        throw new Error(`Failed to solve saddle point system: ${error.message}`);
    }
}

/**
 * Projects gradient onto constraint tangent space (Section 5.3.1)
 * Solves the optimization problem:
 * min_g̃ ||g̃ - g||²_H^s_γ s.t. Cg̃ = 0
 *
 * @param {Array} gradient - Flattened gradient vector
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Object} constraintData - Constraint Jacobian and other data
 * @param {Object} params - Additional parameters (alpha, beta)
 * @returns {Array} - Projected gradient
 */
export function projectGradient(gradient, vertices, edges, constraintData, params) {
    console.log("==== PROJECT GRADIENT ONTO CONSTRAINT TANGENT SPACE ====");
    console.log("Following Section 5.3.1 of the paper");
    
    // Extract constraint Jacobian
    const { jacobian } = constraintData;
    
    if (!jacobian || jacobian.length === 0) {
        console.log("No constraints to project against, returning original gradient");
        console.log("======================================================");
        return gradient;
    }
    
    try {
        // Build the inner product matrix
        console.log("Building inner product matrix A_bar");
        const verticesCopy = vertices.map(v => [...v]);
        const result = build_A_bar_2D(
            params.alpha,
            params.beta,
            verticesCopy,
            edges
        );
        
        // Extract A_bar
        const A_bar = result.A_bar.toArray();
        
        /*
         * For gradient projection, we want to solve:
         * min_g̃ ||g̃ - g||²_H^s_γ s.t. Cg̃ = 0
         * 
         * This leads to the saddle point system:
         * [ A C^T ] [ g̃ ] = [ Ag ]
         * [ C  0  ] [ λ ]   [ 0  ]
         * 
         * Where A is the inner product matrix, C is the constraint Jacobian,
         * g̃ is the projected gradient we're seeking, and λ are the Lagrange multipliers.
         */
        
        // Compute Ag (A * gradient)
        console.log("Computing Ag for the right-hand side");
        const Ag = math.multiply(A_bar, gradient);
        
        // Solve the saddle point system
        console.log("Solving saddle point system for gradient projection");
        const { x: projectedGradient } = solveSaddlePointSystem(
            A_bar,
            jacobian,
            Ag,
            new Array(jacobian.length).fill(0) // d = 0 for gradient projection
        );
        
        // Verify constraint satisfaction
        const constraintViolation = math.norm(math.multiply(jacobian, projectedGradient));
        console.log(`Constraint violation after projection: ${constraintViolation.toExponential(4)}`);
        
        // Compare with original gradient
        const originalNorm = math.norm(gradient);
        const projectedNorm = math.norm(projectedGradient);
        console.log(`Original gradient norm: ${originalNorm.toExponential(4)}`);
        console.log(`Projected gradient norm: ${projectedNorm.toExponential(4)}`);
        console.log(`Relative change: ${((projectedNorm - originalNorm) / originalNorm).toExponential(4)}`);
        
        console.log("======================================================");
        return projectedGradient;
    } catch (error) {
        console.error(`Error in gradient projection: ${error.message}`);
        console.log("======================================================");
        return gradient; // Return original gradient as fallback
    }
}

/**
 * Projects curve back onto constraint set (Section 5.3.2)
 * Solves the optimization problem:
 * min_x ||x||²_H^s_γ s.t. Cx = -Φ(γ̃)
 *
 * @param {Array} vertices - Vertex positions after taking step
 * @param {Array} edges - Edge connections
 * @param {Object} constraints - Constraint definitions
 * @param {Object} constraintData - Constraint values and Jacobian
 * @param {Object} params - Additional parameters (alpha, beta)
 * @returns {Array} - Projected vertices
 */
export function projectOntoConstraintSet(
    vertices, 
    edges,
    constraints, 
    constraintData, 
    params
) {
    console.log("==== PROJECT ONTO CONSTRAINT SET ====");
    console.log("Following Section 5.3.2 of the paper");
    
    if (!constraints || Object.keys(constraints).filter(k => constraints[k]).length === 0) {
        console.log("No active constraints defined, skipping projection");
        console.log("====================================");
        return vertices.map(v => [...v]);
    }
    
    // Extract constraint values and Jacobian
    const { values, evaluate, buildJacobian } = constraintData;
    
    if (!values || values.length === 0) {
        console.log("No constraint values to project, returning original vertices");
        console.log("====================================");
        return vertices.map(v => [...v]);
    }
    
    console.log(`Initial constraint values: [${values.map(v => v.toExponential(4)).join(', ')}]`);
    console.log(`Initial constraint violation: ${math.norm(values).toExponential(4)}`);
    
    // Working copy of vertices
    let projectedVertices = vertices.map(v => [...v]);
    
    // Settings from config
    const tolerance = ProjectionNumerics.CONSTRAINT_TOL;
    const maxIterations = get(config).maxConstraintIterations || 10;
    
    // Track best solution
    let bestVertices = projectedVertices;
    let bestViolation = math.norm(values);
    
    console.log(`Will perform at most ${maxIterations} projection iterations with tolerance ${tolerance}`);
    
    // Newton-like iterations as described in Section 5.3.2
    for (let iter = 0; iter < maxIterations; iter++) {
        console.log(`==== Projection iteration ${iter+1}/${maxIterations} ====`);
        
        // Re-evaluate constraint values at current position
        const constraintValues = evaluate(projectedVertices, constraints);
        const constraintNorm = math.norm(constraintValues);
        
        console.log(`Constraint values: [${constraintValues.map(v => v.toExponential(4)).join(', ')}]`);
        console.log(`Constraint violation: ${constraintNorm.toExponential(4)}`);
        
        // Update best solution if this one is better
        if (constraintNorm < bestViolation) {
            bestViolation = constraintNorm;
            bestVertices = projectedVertices.map(v => [...v]);
            console.log(`New best solution with violation: ${bestViolation.toExponential(4)}`);
        }
        
        // Check if we're close enough
        if (constraintNorm < tolerance) {
            console.log("Constraint projection converged within tolerance");
            break;
        }
        
        // Update Jacobian at current position
        console.log("Updating constraint Jacobian at current position");
        const currentJacobian = buildJacobian(projectedVertices, constraints);
        
        try {
            // Build the inner product matrix for the current position
            console.log("Building inner product matrix for constraint projection");
            const result = build_A_bar_2D(
                params.alpha,
                params.beta,
                projectedVertices,
                edges
            );
            
            const A_bar = result.A_bar.toArray();
            
            /*
             * For constraint projection, we solve the optimization problem:
             * min_x ||x||²_H^s_γ s.t. Cx = -Φ(γ̃)
             * 
             * This leads to the saddle point system:
             * [ A C^T ] [ x ] = [ 0 ]
             * [ C  0  ] [ μ ] = [-Φ(γ̃)]
             * 
             * Where x is the correction to apply to γ̃ to satisfy the constraints.
             */
            
            // Solve the saddle point system
            console.log("Solving saddle point system for constraint projection");
            const { x: correction } = solveSaddlePointSystem(
                A_bar,
                currentJacobian,
                new Array(projectedVertices.length * 2).fill(0), // b = 0
                constraintValues.map(v => -v) // d = -Φ(γ̃)
            );
            
            const correctionNorm = math.norm(correction);
            console.log(`Correction norm: ${correctionNorm.toExponential(4)}`);
            
            // Apply correction with damping if needed
            let dampingFactor = 1.0;
            
            // If correction is too large, apply damping
            if (correctionNorm > ProjectionNumerics.MAX_CORRECTION_NORM) {
                const oldFactor = dampingFactor;
                dampingFactor = ProjectionNumerics.MAX_CORRECTION_NORM / correctionNorm;
                console.log(`Reducing damping factor due to large correction: ${oldFactor} -> ${dampingFactor}`);
            }
            
            // Always apply some damping for stability
            if (ProjectionNumerics.DAMPING_FACTOR < 1.0) {
                dampingFactor *= ProjectionNumerics.DAMPING_FACTOR;
                console.log(`Applied additional damping, final factor: ${dampingFactor}`);
            }
            
            const verticesBeforeCorrection = projectedVertices.map(v => [...v]);
            
            // Apply the correction to vertices
            for (let i = 0; i < projectedVertices.length; i++) {
                projectedVertices[i][0] += correction[i * 2] * dampingFactor;
                projectedVertices[i][1] += correction[i * 2 + 1] * dampingFactor;
            }
            
            // Check if correction made things worse
            const newConstraintValues = evaluate(projectedVertices, constraints);
            const newConstraintNorm = math.norm(newConstraintValues);
            
            console.log(`After correction, constraint violation: ${newConstraintNorm.toExponential(4)}`);
            
            if (newConstraintNorm > constraintNorm * 1.1) { // Allow slight increase to avoid getting stuck
                console.log(`WARNING: Correction increased violation (${constraintNorm.toExponential(4)} -> ${newConstraintNorm.toExponential(4)}), reverting`);
                
                // Try with smaller step
                const smallerDampingFactor = dampingFactor * 0.5;
                console.log(`Trying smaller damping factor: ${dampingFactor} -> ${smallerDampingFactor}`);
                
                // Reset and try smaller step
                projectedVertices = verticesBeforeCorrection.map(v => [...v]);
                
                // Apply the correction with smaller damping
                for (let i = 0; i < projectedVertices.length; i++) {
                    projectedVertices[i][0] += correction[i * 2] * smallerDampingFactor;
                    projectedVertices[i][1] += correction[i * 2 + 1] * smallerDampingFactor;
                }
                
                // Check again
                const newSmallConstraintValues = evaluate(projectedVertices, constraints);
                const newSmallConstraintNorm = math.norm(newSmallConstraintValues);
                
                console.log(`With smaller damping, constraint violation: ${newSmallConstraintNorm.toExponential(4)}`);
                
                if (newSmallConstraintNorm > constraintNorm) {
                    console.log(`Smaller correction still increased violation, reverting to previous state`);
                    projectedVertices = verticesBeforeCorrection;
                    break; // Stop iterations if we can't make progress
                }
            }
            
        } catch (error) {
            console.error(`Error in constraint projection iteration: ${error.message}`);
            projectedVertices = bestVertices;
            break;
        }
    }
    
    // Return best solution
    const finalConstraintValues = evaluate(projectedVertices, constraints);
    const finalViolation = math.norm(finalConstraintValues);
    console.log(`Final constraint values: [${finalConstraintValues.map(v => v.toExponential(4)).join(', ')}]`);
    console.log(`Final constraint violation: ${finalViolation.toExponential(4)}`);
    
    // If best solution is significantly better, use that instead
    if (finalViolation > bestViolation * 1.1) {
        console.log(`Final solution (violation=${finalViolation.toExponential(4)}) is worse than best found (violation=${bestViolation.toExponential(4)}), using best found`);
        projectedVertices = bestVertices;
    }
    
    console.log("====================================");
    return projectedVertices;
}