// src/lib/constraintProjection.js
import * as math from 'mathjs';
import { build_A_bar_2D } from './innerProduct';
import { get } from 'svelte/store';
import { config } from './stores';

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

        // Add regularization to stabilize the system
        // For A_bar, add a small diagonal term
        const epsilon = get(config).epsilonStability || 1e-7;
        const regularizedA = math.add(A_mat, math.multiply(epsilon, math.identity(n)));
        
        // Build the full saddle point matrix
        // [ A_bar   C^T ]
        // [ C      0   ]
        const topRow = math.concat(regularizedA, C_transpose, 1);
        const bottomRow = math.concat(C_mat, math.multiply(-epsilon, math.identity(m)), 1); // Use -ϵI instead of 0
        const fullMatrix = math.concat(topRow, bottomRow, 0);
        
        // Build full right-hand side vector
        const fullRHS = math.concat(b_vec, d_vec);
        
        // Solve the system
        console.log("Solving linear system...");
        const solution = math.lusolve(fullMatrix, fullRHS);
        
        // Extract solution components
        const x = math.subset(solution, math.index(math.range(0, n), 0)).toArray().flat();
        const lambda = math.subset(solution, math.index(math.range(n, n + m), 0)).toArray().flat();
        
        console.log("Solution found: |x| =", math.norm(x), ", |λ| =", math.norm(lambda));
        
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
    console.log("Gradient:", gradient.slice(0, 10) + (gradient.length > 10 ? "..." : ""));
    console.log("Vertices length:", vertices.length);
    console.log("Raw edges:", edges);
    
    // CRITICAL FIX: Ensure edges are properly structured before proceeding
    if (!Array.isArray(edges)) {
        console.error("ERROR: edges is not an array:", edges);
        console.log("===== PROJECT GRADIENT - END WITH ERROR =====");
        throw new Error(`Invalid edges parameter (not an array): ${typeof edges}`);
    }
    
    // Make a deep copy of the edges array to prevent mutations
    const safeEdges = edges.map(e => Array.isArray(e) ? [...e] : e);
    console.log("Safe copy of edges:", JSON.stringify(safeEdges));
    
    // Extract constraint Jacobian and other needed data
    const { jacobian } = constraintData;
    
    if (!jacobian || jacobian.length === 0) {
        console.log("No constraints to project against, returning original gradient");
        console.log("===== PROJECT GRADIENT - END =====");
        return gradient;
    }
    
    // Build the fractional Sobolev inner product matrix (A_bar)
    console.log("Building A_bar matrix for gradient projection");
    console.log("Calling build_A_bar_2D with:", 
        `alpha=${params.alpha}`,
        `beta=${params.beta}`,
        `vertices (length ${vertices.length})`,
        `safeEdges (length ${safeEdges.length})`
    );
    
    try {
        // Pass safe copies of all parameters
        const verticesCopy = vertices.map(v => [...v]);
        
        // CRITICAL FIX: Match parameter order with function definition
        // From innerProduct.js: build_A_bar_2D(alpha, beta, vertices, edges)
        const result = build_A_bar_2D(
            params.alpha,
            params.beta,
            verticesCopy,
            safeEdges
        );
        
        // Extract A_bar from the result object
        const A_bar = result.A_bar;
        
        console.log("A_bar matrix built successfully:", A_bar.size()[0], "×", A_bar.size()[1]);
        
        // Convert A_bar matrix to array for compatibility with solveSaddlePointSystem
        const A_bar_array = A_bar.toArray();
        
        // Solve the saddle point system
        console.log("Solving saddle point system for gradient projection");
        const { x: projectedGradient } = solveSaddlePointSystem(
            A_bar_array,
            jacobian,
            gradient,
            new Array(jacobian.length).fill(0) // d = 0 for gradient projection
        );
        
        // CRITICAL FIX: Check for numerical instability in the result
        const gradientNorm = math.norm(projectedGradient);
        console.log("Projected gradient norm:", gradientNorm);
        
        // Apply clamping if the gradient is too large
        const MAX_GRADIENT_NORM = 10.0;  // Adjust this value based on your application scale
        let normalizedGradient = projectedGradient;
        
        if (gradientNorm > MAX_GRADIENT_NORM) {
            console.log(`Normalizing large gradient (${gradientNorm} -> ${MAX_GRADIENT_NORM})`);
            const scaleFactor = MAX_GRADIENT_NORM / gradientNorm;
            normalizedGradient = projectedGradient.map(val => val * scaleFactor);
        }
        
        console.log("Final projected gradient norm:", math.norm(normalizedGradient));
        console.log("===== PROJECT GRADIENT - END =====");
        return normalizedGradient;
    } catch (error) {
        console.error("Error building A_bar or solving system:", error);
        console.log("===== PROJECT GRADIENT - END WITH ERROR =====");
        throw error;
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
    console.log("Projecting curve onto constraint set");
    console.log("Vertices length:", vertices.length);
    console.log("Edges:", JSON.stringify(edges));
    console.log("Alpha:", params.alpha, "Beta:", params.beta);
    
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
    
    console.log("Initial constraint violation:", math.norm(values));
    
    // Make a working copy of vertices
    let projectedVertices = vertices.map(v => [...v]);
    
    const tolerance = get(config).constraintTolerance || 1e-4;
    const maxIterations = get(config).maxConstraintIterations || 3;
    
    // CRITICAL FIX: Add damping to prevent too large corrections
    let dampingFactor = 0.3; // Start with a conservative value
    
    // Newton-like iterations to project onto constraint set
    for (let iter = 0; iter < maxIterations; iter++) {
        // Re-evaluate constraint values at current point
        const constraintValues = evaluate(projectedVertices, constraints);
        const constraintNorm = math.norm(constraintValues);
        
        console.log(`Iteration ${iter+1}/${maxIterations}, constraint violation: ${constraintNorm}`);
        
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
            // CRITICAL FIX: Match parameter order with function definition
            // From innerProduct.js: build_A_bar_2D(alpha, beta, vertices, edges)
            const result = build_A_bar_2D(
                params.alpha,
                params.beta,
                projectedVertices, 
                safeEdges
            );
            
            // Extract A_bar from the result object
            const A_bar = result.A_bar;
            
            console.log("A_bar matrix built successfully:", A_bar.size()[0], "×", A_bar.size()[1]);
            
            // Convert A_bar matrix to array for compatibility with solveSaddlePointSystem
            const A_bar_array = A_bar.toArray();
            
            // Solve the saddle point system to find the correction
            console.log("Solving saddle point system for constraint projection");
            const { x: correction } = solveSaddlePointSystem(
                A_bar_array,
                currentJacobian,
                new Array(projectedVertices.length * 2).fill(0), // b = 0
                constraintValues.map(v => -v) // d = -Φ(γ̃)
            );
            
            // Apply correction with damping and clamping
            const correctionNorm = math.norm(correction);
            console.log(`Raw correction magnitude: ${correctionNorm}`);
            
            // CRITICAL FIX: Apply maximum correction limit
            const MAX_CORRECTION = 100.0; // Adjust based on your scene scale
            let scaledCorrection = correction;
            
            if (correctionNorm > MAX_CORRECTION) {
                const scaleFactor = MAX_CORRECTION / correctionNorm;
                scaledCorrection = correction.map(val => val * scaleFactor);
                console.log(`Limiting large correction: ${correctionNorm} -> ${MAX_CORRECTION}`);
            }
            
            // Apply damping
            const dampedCorrection = scaledCorrection.map(val => val * dampingFactor);
            const dampedNorm = math.norm(dampedCorrection);
            console.log(`Applying damped correction with magnitude: ${dampedNorm} (damping=${dampingFactor})`);
            
            // Apply the correction to vertices
            const verticesBeforeCorrection = projectedVertices.map(v => [...v]); // Save for backtracking
            
            // Reshape correction back to vertex format
            for (let i = 0; i < projectedVertices.length; i++) {
                projectedVertices[i][0] += dampedCorrection[i * 2];
                projectedVertices[i][1] += dampedCorrection[i * 2 + 1];
            }
            
            // CRITICAL FIX: Backtracking if correction made things worse
            const newConstraintValues = evaluate(projectedVertices, constraints);
            const newConstraintNorm = math.norm(newConstraintValues);
            
            if (newConstraintNorm > constraintNorm) {
                console.log(`Correction increased violation (${constraintNorm} -> ${newConstraintNorm}), backtracking`);
                projectedVertices = verticesBeforeCorrection;
                
                // Reduce damping for next iteration to be more conservative
                if (iter < maxIterations - 1) {
                    dampingFactor *= 0.5;
                }
            }
            
        } catch (error) {
            console.error("Error building A_bar or solving system:", error);
            console.log("===== PROJECT ONTO CONSTRAINT SET - END WITH ERROR =====");
            throw error;
        }
    }
    
    // Final evaluation of constraints
    const finalViolation = math.norm(evaluate(projectedVertices, constraints));
    console.log(`Final constraint violation after projection: ${finalViolation}`);
    console.log("===== PROJECT ONTO CONSTRAINT SET - END =====");
    
    return projectedVertices;
}