// src/lib/constraintProjection.js
import * as math from 'mathjs';
import { build_A_bar_2D } from './innerProduct';
import { get } from 'svelte/store';
import { config } from './stores';

// Configuration options for numerical stability measures
// Set these to false to follow the paper exactly without stabilization
const STABILIZATION_CONFIG = {
    // Whether to apply any stabilization measures (master switch)
    useStabilization: false,
    
    // Saddle point system regularization
    useRegularization: true,      // Add small diagonal terms to stabilize matrix inversion
    regularizationEpsilon: 1e-6,  // Regularization parameter (default)
    
    // Gradient projection stabilization
    limitGradientMagnitude: true, // Cap gradient magnitude to prevent extreme steps
    maxGradientNorm: 1000.0,      // Maximum allowed gradient norm if limiting is enabled
    
    // Constraint projection stabilization 
    useDamping: true,             // Apply damping factor to corrections
    dampingFactor: 0.3,           // Initial damping factor for corrections (0-1)
    limitCorrectionMagnitude: true, // Cap correction magnitude
    maxCorrectionNorm: 100.0,      // Maximum allowed correction if limiting is enabled
    useBacktracking: true,        // Discard corrections that increase constraint violation
};

/**
 * Solves the saddle point system for constrained optimization
 * 
 * [ A_bar   C^T ] [ x ] = [ b ]
 * [ C      0   ] [ λ ] = [ d ]
 *
 * where:
 * - A_bar is the (2|V| × 2|V|) fractional Sobolev inner product matrix
 * - C is the Jacobian of constraints (|C| × 2|V|)
 * - x is the solution vector in R^(2|V|)
 * - λ are the Lagrange multipliers
 * - b and d are the right-hand side vectors
 *
 * @param {Array|Matrix} A_bar - Fractional Sobolev inner product matrix
 * @param {Array|Matrix} C - Constraint Jacobian matrix
 * @param {Array} b - First right-hand side vector
 * @param {Array} d - Second right-hand side vector
 * @returns {Object} - Solution vectors x and λ
 */
export function solveSaddlePointSystem(A_bar, C, b, d) {
    console.log("Solving saddle point system");
    console.log("A_bar dimensions:", A_bar.length, "×", A_bar[0]?.length || 0);
    console.log("C dimensions:", C.length, "×", C[0]?.length || 0);
    console.log("b length:", b.length);
    console.log("d length:", d.length);

    try {
        // Convert arrays to math.js matrices for better handling
        const A_mat = math.matrix(A_bar);
        const C_mat = math.matrix(C);
        const C_transpose = math.transpose(C_mat);
        const b_vec = math.matrix(b);
        const d_vec = math.matrix(d);

        // Get dimensions
        const n = b.length; // 2|V|
        const m = d.length; // Number of constraints
        
        console.log(`System dimensions: n=${n}, m=${m}`);

        // Build the saddle point matrix with regularization
        let regularizedA = A_mat;
        let zeroBlock = math.zeros(m, m);
        
        // Apply regularization to improve conditioning
        if (STABILIZATION_CONFIG.useStabilization && STABILIZATION_CONFIG.useRegularization) {
            // Get configured epsilon or use default
            const epsilon = get(config).epsilonStability || 
                          STABILIZATION_CONFIG.regularizationEpsilon;
                          
            console.log(`Applying regularization with epsilon=${epsilon}`);
            
            // Add εI to A_bar block
            regularizedA = math.add(A_mat, math.multiply(epsilon, math.identity(n)));
            
            // Use -εI instead of 0 for the bottom-right block
            zeroBlock = math.multiply(-epsilon, math.identity(m));
        }
        
        // [ A_bar   C^T ]
        // [ C      0   ]
        const topRow = math.concat(regularizedA, C_transpose, 1);
        const bottomRow = math.concat(C_mat, zeroBlock, 1);
        const fullMatrix = math.concat(topRow, bottomRow, 0);
        
        // Build full right-hand side vector
        const fullRHS = math.concat(b_vec, d_vec);
        
        // Solve the system
        console.log("Solving linear system...");
        const solution = math.lusolve(fullMatrix, fullRHS);
        
        // Extract solution components
        const x = math.subset(solution, math.index(math.range(0, n), 0)).toArray().flat();
        const lambda = math.subset(solution, math.index(math.range(n, n + m), 0)).toArray().flat();
        
        const xNorm = math.norm(x);
        const lambdaNorm = math.norm(lambda);
        
        console.log(`Solution found: |x| = ${xNorm}, |λ| = ${lambdaNorm}`);
        
        // Check for extremely large solutions that indicate numerical issues
        if (!isFinite(xNorm) || xNorm > 1e10) {
            console.error("WARNING: Solution has extremely large magnitude, likely numerical instability");
        }
        
        return { x, lambda };
    } catch (error) {
        console.error("Error solving saddle point system:", error);
        throw error; // Rethrow the error - never silently recover
    }
}

/**
 * Projects gradient onto constraint tangent space (Section 5.3.1)
 * Solving the saddle point system:
 * 
 * [ A_bar   C^T ] [ g̃ ] = [ dE^α_β|ᵧ^T ]
 * [ C       0  ] [ λ ]   [     0      ]
 *
 * @param {Array} gradient - Flattened gradient vector (2|V| x 1)
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Object} constraintData - Constraint Jacobian and other data
 * @param {Object} params - Additional parameters (alpha, beta)
 * @returns {Array} - Projected gradient (2|V| x 1)
 */
export function projectGradient(gradient, vertices, edges, constraintData, params) {
    console.log("===== PROJECT GRADIENT - START =====");
    console.log("Projecting gradient using full saddle point system");
    
    // Log gradient norm for diagnosis
    const gradientNorm = math.norm(gradient);
    console.log(`Original gradient norm: ${gradientNorm}`);
    
    // Ensure edges are properly structured before proceeding
    if (!Array.isArray(edges)) {
        console.error("ERROR: edges is not an array:", edges);
        console.log("===== PROJECT GRADIENT - END WITH ERROR =====");
        throw new Error(`Invalid edges parameter (not an array): ${typeof edges}`);
    }
    
    // Make a deep copy of the edges array to prevent mutations
    const safeEdges = edges.map(e => Array.isArray(e) ? [...e] : e);
    
    // Extract constraint Jacobian and other needed data
    const { jacobian } = constraintData;
    
    if (!jacobian || jacobian.length === 0) {
        console.log("No constraints to project against, returning original gradient");
        console.log("===== PROJECT GRADIENT - END =====");
        return gradient;
    }
    
    try {
        // Build the fractional Sobolev inner product matrix (A_bar)
        console.log("Building A_bar matrix for gradient projection");
        
        // Pass safe copies of all parameters
        const verticesCopy = vertices.map(v => [...v]);
        
        // Build A_bar matrix
        const result = build_A_bar_2D(
            params.alpha,
            params.beta,
            verticesCopy,
            safeEdges
        );
        
        // Extract A_bar from the result object
        const A_bar = result.A_bar;
        const A_bar_array = A_bar.toArray();
        
        // Solve the saddle point system
        console.log("Solving saddle point system for gradient projection");
        const { x: projectedGradient } = solveSaddlePointSystem(
            A_bar_array,
            jacobian,
            gradient,
            new Array(jacobian.length).fill(0) // d = 0 for gradient projection
        );
        
        // Check for numerical instability in the result and limit if configured
        let finalGradient = projectedGradient;
        const projGradNorm = math.norm(projectedGradient);
        
        console.log(`Projected gradient norm: ${projGradNorm}`);
        
        // Apply clamping if enabled and the gradient is too large
        if (STABILIZATION_CONFIG.useStabilization && 
            STABILIZATION_CONFIG.limitGradientMagnitude && 
            projGradNorm > STABILIZATION_CONFIG.maxGradientNorm) {
            
            console.log(`Normalizing large gradient (${projGradNorm} -> ${STABILIZATION_CONFIG.maxGradientNorm})`);
            const scaleFactor = STABILIZATION_CONFIG.maxGradientNorm / projGradNorm;
            finalGradient = projectedGradient.map(val => val * scaleFactor);
            console.log(`Final projected gradient norm after scaling: ${math.norm(finalGradient)}`);
        }
        
        console.log("===== PROJECT GRADIENT - END =====");
        return finalGradient;
    } catch (error) {
        console.error("Error in gradient projection:", error);
        console.log("===== PROJECT GRADIENT - END WITH ERROR =====");
        // Fallback to original gradient if projection fails
        return gradient;
    }
}

/**
 * Projects curve back onto constraint set (Section 5.3.2)
 * Solves: min_x (1/2)x^T A̅ x s.t. Cx = -Φ(γ̃)
 * Through the saddle point system:
 * 
 * [ A_bar   C^T ] [ x ] = [    0     ]
 * [ C       0  ] [ μ ]   [ -Φ(γ̃)    ]
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

    // Verify edges is properly defined and is an array
    if (!Array.isArray(edges)) {
        console.error("ERROR: edges is not an array:", edges);
        console.log("===== PROJECT ONTO CONSTRAINT SET - END WITH ERROR =====");
        throw new Error(`Invalid edges parameter (not an array): ${typeof edges}`);
    }
    
    // Make a deep copy of the edges array to prevent mutations
    const safeEdges = edges.map(e => Array.isArray(e) ? [...e] : e);
    
    const initialViolation = math.norm(values);
    console.log("Initial constraint violation:", initialViolation);
    
    // Make a working copy of vertices
    let projectedVertices = vertices.map(v => [...v]);
    
    const tolerance = get(config).constraintTolerance || 1e-4;
    const maxIterations = get(config).maxConstraintIterations || 3;
    
    // Set damping factor if enabled
    let dampingFactor = STABILIZATION_CONFIG.useStabilization && STABILIZATION_CONFIG.useDamping ? 
                        STABILIZATION_CONFIG.dampingFactor : 1.0;
    
    // Track best solution so far
    let bestVertices = projectedVertices;
    let bestViolation = initialViolation;
    
    // Newton-like iterations to project onto constraint set
    for (let iter = 0; iter < maxIterations; iter++) {
        // Re-evaluate constraint values at current point
        const constraintValues = evaluate(projectedVertices, constraints);
        const constraintNorm = math.norm(constraintValues);
        
        console.log(`Iteration ${iter+1}/${maxIterations}, constraint violation: ${constraintNorm}`);
        
        // Update best solution if this one is better
        if (constraintNorm < bestViolation) {
            bestViolation = constraintNorm;
            bestVertices = projectedVertices.map(v => [...v]);
        }
        
        // Check if we're already close enough to the constraint set
        if (constraintNorm < tolerance) {
            console.log("Constraint projection converged");
            break;
        }
        
        // Update Jacobian at current position
        const currentJacobian = buildJacobian(projectedVertices, constraints);
        
        // Build the fractional Sobolev inner product matrix (A_bar)
        console.log("Building A_bar matrix for constraint projection");
        
        try {
            // Match parameter order with function definition
            const result = build_A_bar_2D(
                params.alpha,
                params.beta,
                projectedVertices, 
                safeEdges
            );
            
            // Extract A_bar from the result object
            const A_bar = result.A_bar;
            const A_bar_array = A_bar.toArray();
            
            // Solve the saddle point system to find the correction
            console.log("Solving saddle point system for constraint projection");
            const { x: correction } = solveSaddlePointSystem(
                A_bar_array,
                currentJacobian,
                new Array(projectedVertices.length * 2).fill(0), // b = 0
                constraintValues.map(v => -v) // d = -Φ(γ̃)
            );
            
            // Apply correction with damping and clamping if enabled
            const correctionNorm = math.norm(correction);
            console.log(`Raw correction magnitude: ${correctionNorm}`);
            
            let finalCorrection = correction;
            
            // Apply maximum correction limit if configured
            if (STABILIZATION_CONFIG.useStabilization && 
                STABILIZATION_CONFIG.limitCorrectionMagnitude && 
                correctionNorm > STABILIZATION_CONFIG.maxCorrectionNorm) {
                
                console.log(`Limiting large correction (${correctionNorm} -> ${STABILIZATION_CONFIG.maxCorrectionNorm})`);
                const scaleFactor = STABILIZATION_CONFIG.maxCorrectionNorm / correctionNorm;
                finalCorrection = correction.map(val => val * scaleFactor);
            }
            
            // Apply damping if enabled
            if (STABILIZATION_CONFIG.useStabilization && STABILIZATION_CONFIG.useDamping) {
                console.log(`Applying damping factor: ${dampingFactor}`);
                finalCorrection = finalCorrection.map(val => val * dampingFactor);
            }
            
            // Save vertices for potential backtracking
            const verticesBeforeCorrection = projectedVertices.map(v => [...v]);
            
            // Apply the correction to vertices
            for (let i = 0; i < projectedVertices.length; i++) {
                projectedVertices[i][0] += finalCorrection[i * 2];
                projectedVertices[i][1] += finalCorrection[i * 2 + 1];
            }
            
            // Apply backtracking if correction made things worse and backtracking is enabled
            if (STABILIZATION_CONFIG.useStabilization && STABILIZATION_CONFIG.useBacktracking) {
                const newConstraintValues = evaluate(projectedVertices, constraints);
                const newConstraintNorm = math.norm(newConstraintValues);
                
                if (newConstraintNorm > constraintNorm) {
                    console.log(`Correction increased violation (${constraintNorm} -> ${newConstraintNorm}), backtracking`);
                    projectedVertices = verticesBeforeCorrection;
                    
                    // Reduce damping for next iteration to be more conservative
                    dampingFactor *= 0.5;
                    console.log(`Reduced damping factor to ${dampingFactor}`);
                }
            }
            
        } catch (error) {
            console.error("Error in constraint projection iteration:", error);
            // Fall back to best solution so far
            projectedVertices = bestVertices;
            break;
        }
    }
    
    // Return best solution found if final iteration wasn't better
    const finalViolation = math.norm(evaluate(projectedVertices, constraints));
    if (finalViolation > bestViolation) {
        console.log(`Final solution (violation=${finalViolation}) is worse than best found (violation=${bestViolation}), returning best solution`);
        projectedVertices = bestVertices;
    }
    
    console.log(`Final constraint violation after projection: ${math.norm(evaluate(projectedVertices, constraints))}`);
    console.log("===== PROJECT ONTO CONSTRAINT SET - END =====");
    
    return projectedVertices;
}