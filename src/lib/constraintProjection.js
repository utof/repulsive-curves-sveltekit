/**
 * Direct implementation of constraint projection without math.js
 * Follows the mathematical formulation in the paper's Section 5.3.1 and 5.3.2
 */

import { build_A_bar_2D } from './innerProduct';
import { get } from 'svelte/store';
import { config, initialTotalLength } from './stores'; // Fixed import to include initialTotalLength directly

/**
 * Solves a linear system Ax = b using Gaussian elimination
 * Simpler but more robust than relying on math.js
 * 
 * @param {Array} A - 2D array representing matrix A
 * @param {Array} b - 1D array representing vector b
 * @returns {Array} - Solution vector x
 */
function solveLinearSystem(A, b) {
    const n = A.length;
    if (n === 0) return [];
    
    // Make copies to avoid modifying the originals
    const augmented = new Array(n);
    for (let i = 0; i < n; i++) {
        augmented[i] = [...A[i], b[i]];
    }
    
    // Forward elimination
    for (let i = 0; i < n; i++) {
        // Find pivot
        let maxRow = i;
        let maxVal = Math.abs(augmented[i][i]);
        
        for (let j = i + 1; j < n; j++) {
            const absVal = Math.abs(augmented[j][i]);
            if (absVal > maxVal) {
                maxVal = absVal;
                maxRow = j;
            }
        }
        
        // Swap rows if needed
        if (maxRow !== i) {
            [augmented[i], augmented[maxRow]] = [augmented[maxRow], augmented[i]];
        }
        
        // Skip if the matrix is singular
        if (Math.abs(augmented[i][i]) < 1e-10) {
            console.warn(`Near-singular matrix detected at row ${i}`);
            augmented[i][i] = 1e-10; // Add small regularization
        }
        
        // Eliminate below
        for (let j = i + 1; j < n; j++) {
            const factor = augmented[j][i] / augmented[i][i];
            for (let k = i; k <= n; k++) {
                augmented[j][k] -= factor * augmented[i][k];
            }
        }
    }
    
    // Back substitution
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
        let sum = 0;
        for (let j = i + 1; j < n; j++) {
            sum += augmented[i][j] * x[j];
        }
        x[i] = (augmented[i][n] - sum) / augmented[i][i];
    }
    
    return x;
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
    console.log("Solving saddle point system");
    console.log("A dimensions:", A.length, "×", A[0]?.length || 0);
    console.log("C dimensions:", C.length, "×", C[0]?.length || 0);
    console.log("b length:", b.length);
    console.log("d length:", d.length);
    
    try {
        const n = A.length;
        const m = C.length;
        
        // For extremely small constraint systems, use a direct method
        if (m === 1) {
            return solveSimplifiedConstraintSystem(A, C[0], b, d[0]);
        }
        
        // Build the saddle point matrix [A C^T; C 0]
        const fullSize = n + m;
        const fullMatrix = Array(fullSize).fill().map(() => Array(fullSize).fill(0));
        
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
        
        // Fill zeros block (possibly with small negative diagonal for stability)
        const epsilon = 1e-10; // Small regularization
        for (let i = 0; i < m; i++) {
            fullMatrix[n + i][n + i] = -epsilon;
        }
        
        // Build full right-hand side vector
        const fullRHS = [...b, ...d];
        
        // Solve the system
        console.log("Solving linear system...");
        const solution = solveLinearSystem(fullMatrix, fullRHS);
        
        // Extract components
        const x = solution.slice(0, n);
        const lambda = solution.slice(n, n + m);
        
        const xNorm = Math.sqrt(x.reduce((sum, val) => sum + val*val, 0));
        const lambdaNorm = Math.sqrt(lambda.reduce((sum, val) => sum + val*val, 0));
        
        console.log("Solution found: |x| =", xNorm, ", |λ| =", lambdaNorm);
        return { x, lambda };
    } catch (error) {
        console.error("Error solving saddle point system:", error);
        throw error;
    }
}

/**
 * For a single constraint (m=1), use a simplified and more stable approach
 * This avoids the full saddle point system which might be ill-conditioned
 */
function solveSimplifiedConstraintSystem(A, c, b, d) {
    console.log("Using simplified approach for single constraint");
    const n = A.length;
    
    // Compute A^(-1)b and A^(-1)c
    // First solve Ax = b
    const A_inv_b = solveLinearSystem(A, b);
    
    // Then solve Ay = c for each row of c
    const A_inv_c = solveLinearSystem(A, c);
    
    // Compute λ = (c^T A^(-1) c)^(-1) (c^T A^(-1) b - d)
    let cTAinvc = 0;
    let cTAinvb = 0;
    
    for (let i = 0; i < n; i++) {
        cTAinvc += c[i] * A_inv_c[i];
        cTAinvb += c[i] * A_inv_b[i];
    }
    
    if (Math.abs(cTAinvc) < 1e-10) {
        console.warn("Near-singular constraint system, adding regularization");
        cTAinvc = Math.sign(cTAinvc) * 1e-10;
    }
    
    const lambda = (cTAinvb - d) / cTAinvc;
    
    // Compute x = A^(-1)b - λ * A^(-1)c
    const x = A_inv_b.map((val, i) => val - lambda * A_inv_c[i]);
    
    return { x, lambda: [lambda] };
}

/**
 * Projects gradient onto constraint tangent space (Section 5.3.1)
 *
 * @param {Array} gradient - Flattened gradient vector
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Object} constraintData - Constraint Jacobian and other data
 * @param {Object} params - Additional parameters (alpha, beta)
 * @returns {Array} - Projected gradient
 */
export function projectGradient(gradient, vertices, edges, constraintData, params) {
    console.log("===== PROJECT GRADIENT - START =====");
    console.log("Projecting gradient onto constraint tangent space");
    
    // Ensure edges are properly structured
    if (!Array.isArray(edges)) {
        console.error("ERROR: edges is not an array:", edges);
        console.log("===== PROJECT GRADIENT - END WITH ERROR =====");
        return gradient; // Return original gradient as fallback
    }
    
    // Make a deep copy of the edges array
    const safeEdges = edges.map(e => Array.isArray(e) ? [...e] : e);
    
    // Extract constraint Jacobian
    const { jacobian } = constraintData;
    
    if (!jacobian || jacobian.length === 0) {
        console.log("No constraints to project against, returning original gradient");
        console.log("===== PROJECT GRADIENT - END =====");
        return gradient;
    }
    
    try {
        // Build the inner product matrix
        console.log("Building inner product matrix for gradient projection");
        const verticesCopy = vertices.map(v => [...v]);
        
        const result = build_A_bar_2D(
            params.alpha,
            params.beta,
            verticesCopy,
            safeEdges
        );
        
        // Extract A_bar and convert to array
        const A_bar = result.A_bar.toArray();
        
        // We want to solve: min_g̃ ||g̃ - g||²_A such that Cg̃ = 0
        // This is equivalent to solving the saddle point system:
        // [ A C^T ] [ g̃ ] = [ Ag ]
        // [ C  0  ] [ λ ]   [ 0  ]
        
        // Compute Ag (A * gradient)
        const Ag = new Array(gradient.length).fill(0);
        for (let i = 0; i < A_bar.length; i++) {
            for (let j = 0; j < gradient.length; j++) {
                Ag[i] += A_bar[i][j] * gradient[j];
            }
        }
        
        // Solve the saddle point system
        console.log("Solving saddle point system for gradient projection");
        const { x: projectedGradient } = solveSaddlePointSystem(
            A_bar,
            jacobian,
            Ag,
            new Array(jacobian.length).fill(0) // Right-hand side is zero
        );
        
        console.log("===== PROJECT GRADIENT - END =====");
        return projectedGradient;
    } catch (error) {
        console.error("Error in gradient projection:", error);
        console.log("===== PROJECT GRADIENT - END WITH ERROR =====");
        return gradient; // Return original gradient as fallback
    }
}

/**
 * Projects curve back onto constraint set (Section 5.3.2)
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
    console.log("===== PROJECT ONTO CONSTRAINT SET - START =====");
    
    if (!constraints || Object.keys(constraints).length === 0) {
        console.log("No constraints defined, skipping projection");
        console.log("===== PROJECT ONTO CONSTRAINT SET - END =====");
        return vertices.map(v => [...v]);
    }
    
    // Extract constraint values and Jacobian
    const { values, jacobian, evaluate, buildJacobian } = constraintData;
    
    if (!values || values.length === 0) {
        console.log("No constraint values to project, returning original vertices");
        console.log("===== PROJECT ONTO CONSTRAINT SET - END =====");
        return vertices.map(v => [...v]);
    }

    // Special handling for length constraint - direct scaling method
    if (constraints.length && !constraints.barycenter && jacobian.length === 1) {
        console.log("Using direct scaling for length constraint");
        return projectOntoLengthConstraint(vertices, edges, constraints);
    }
    
    // Make a deep copy of the edges array
    const safeEdges = edges.map(e => Array.isArray(e) ? [...e] : e);
    
    // Working copy of vertices
    let projectedVertices = vertices.map(v => [...v]);
    
    // Settings
    const tolerance = get(config).constraintTolerance || 1e-4;
    const maxIterations = get(config).maxConstraintIterations || 3;
    
    // Track best solution
    let bestVertices = projectedVertices;
    let bestViolation = Math.sqrt(values.reduce((sum, v) => sum + v*v, 0));
    
    console.log("Initial constraint violation:", bestViolation);
    
    // Newton-like iterations
    for (let iter = 0; iter < maxIterations; iter++) {
        // Re-evaluate constraint values
        const constraintValues = evaluate(projectedVertices, constraints);
        const constraintNorm = Math.sqrt(constraintValues.reduce((sum, v) => sum + v*v, 0));
        
        console.log(`Iteration ${iter+1}/${maxIterations}, constraint violation: ${constraintNorm}`);
        
        // Update best solution if this one is better
        if (constraintNorm < bestViolation) {
            bestViolation = constraintNorm;
            bestVertices = projectedVertices.map(v => [...v]);
        }
        
        // Check if we're close enough
        if (constraintNorm < tolerance) {
            console.log("Constraint projection converged");
            break;
        }
        
        // Update Jacobian at current position
        const currentJacobian = buildJacobian(projectedVertices, constraints);
        
        try {
            // Build the inner product matrix
            console.log("Building inner product matrix for constraint projection");
            const result = build_A_bar_2D(
                params.alpha,
                params.beta,
                projectedVertices, 
                safeEdges
            );
            
            const A_bar = result.A_bar.toArray();
            
            // Solve the saddle point system
            console.log("Solving saddle point system for constraint projection");
            const { x: correction } = solveSaddlePointSystem(
                A_bar,
                currentJacobian,
                new Array(projectedVertices.length * 2).fill(0), // b = 0
                constraintValues.map(v => -v) // d = -Φ(γ̃)
            );
            
            console.log(`Raw correction magnitude: ${Math.sqrt(correction.reduce((sum, v) => sum + v*v, 0))}`);
            
            // Apply correction (with damping if needed)
            const dampingFactor = 0.5; // Conservative damping factor
            
            const verticesBeforeCorrection = projectedVertices.map(v => [...v]);
            
            // Apply the correction to vertices
            for (let i = 0; i < projectedVertices.length; i++) {
                projectedVertices[i][0] += correction[i * 2] * dampingFactor;
                projectedVertices[i][1] += correction[i * 2 + 1] * dampingFactor;
            }
            
            // Check if correction made things worse
            const newConstraintValues = evaluate(projectedVertices, constraints);
            const newConstraintNorm = Math.sqrt(newConstraintValues.reduce((sum, v) => sum + v*v, 0));
            
            if (newConstraintNorm > constraintNorm) {
                console.log(`Correction increased violation (${constraintNorm} -> ${newConstraintNorm}), backtracking`);
                projectedVertices = verticesBeforeCorrection;
                break; // Stop iterations if we can't make progress
            }
            
        } catch (error) {
            console.error("Error in constraint projection iteration:", error);
            projectedVertices = bestVertices;
            break;
        }
    }
    
    // Return best solution
    const finalConstraintValues = evaluate(projectedVertices, constraints);
    const finalViolation = Math.sqrt(finalConstraintValues.reduce((sum, v) => sum + v*v, 0));
    console.log(`Final constraint violation after projection: ${finalViolation}`);
    
    // If best solution is significantly better, use that instead
    if (finalViolation > bestViolation * 1.1) {
        console.log(`Final solution (violation=${finalViolation}) is worse than best found (violation=${bestViolation}), returning best solution`);
        projectedVertices = bestVertices;
    }
    
    console.log("===== PROJECT ONTO CONSTRAINT SET - END =====");
    return projectedVertices;
}

/**
 * Projects vertices onto length constraint using direct scaling
 * This is much more stable than using the general saddle point system for length constraints
 * 
 * @param {Array} vertices - Vertex positions 
 * @param {Array} edges - Edge connections
 * @param {Object} constraints - Constraint definitions
 * @returns {Array} - Projected vertices
 */
function projectOntoLengthConstraint(vertices, edges, constraints) {
    console.log("Projecting directly onto length constraint");
    
    // Calculate current total length
    let currentLength = 0;
    for (const [v1, v2] of edges) {
        const dx = vertices[v2][0] - vertices[v1][0];
        const dy = vertices[v2][1] - vertices[v1][1];
        currentLength += Math.sqrt(dx*dx + dy*dy);
    }
    
    // Determine target length
    let targetLength = constraints.lengthTarget || 0;
    if (constraints.lengthPercentage) {
        // Fix: Get initialTotalLength directly from the store
        const initialLengthValue = get(initialTotalLength); 
        targetLength = initialLengthValue * (constraints.lengthPercentage / 100);
        console.log(`Using ${constraints.lengthPercentage}% of initial length (${initialLengthValue}) as target: ${targetLength}`);
    }
    
    console.log(`Current length: ${currentLength}, target length: ${targetLength}`);
    
    // If already very close, no need to project
    if (Math.abs(currentLength - targetLength) < 1e-6) {
        console.log("Already at target length, no projection needed");
        return vertices.map(v => [...v]);
    }
    
    // Calculate barycenter (for scaling about the center)
    const barycenter = [0, 0];
    for (const vertex of vertices) {
        barycenter[0] += vertex[0] / vertices.length;
        barycenter[1] += vertex[1] / vertices.length;
    }
    
    // Scale factor
    const scaleFactor = targetLength / currentLength;
    console.log(`Scaling by factor: ${scaleFactor}`);
    
    // Apply scaling about the barycenter
    return vertices.map(v => [
        barycenter[0] + (v[0] - barycenter[0]) * scaleFactor,
        barycenter[1] + (v[1] - barycenter[1]) * scaleFactor
    ]);
}