// src/lib/constraintProjection.js
/**
 * Implements constraint projection using the saddle point system
 * Following Sections 5.3.1 and 5.3.2 of the paper
 * Generalized for both 2D and 3D
 */

import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from './stores';
import { build_A_bar } from './innerProduct';

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
         */
        
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
        
        // Build full right-hand side vector
        const fullRHS = b.concat(d);
        
        // Solve the system
        console.log("Solving saddle point system using math.js...");
        
        const solution = math.lusolve(fullMatrix, fullRHS);
        
        // Extract x and lambda components
        const solutionArray = math.matrix(solution).toArray();
        
        // Extract x and lambda components
        const x = solutionArray.slice(0, n).map(row => row[0]); // Get first column of each row
        const lambda = solutionArray.slice(n, n + m).map(row => row[0]);
        
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
 * Works for both 2D and 3D
 *
 * @param {Array} gradient - Flattened gradient vector
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Object} constraintData - Constraint Jacobian and other data
 * @param {Object} params - Additional parameters (alpha, beta)
 * @returns {Array} - Projected gradient
 */
export function projectGradient(gradient, vertices, edges, constraintData, params) {
    const dimension = get(config).dimension;
    console.log("==== PROJECT GRADIENT ONTO CONSTRAINT TANGENT SPACE ====");
    console.log(`Current dimension: ${dimension}, vertices: ${vertices.length}`);
    console.log(`Gradient array length: ${gradient.length} (expected: ${vertices.length * dimension})`);
    console.log("Following Section 5.3.1 of the paper");
    
    // Extract constraint Jacobian
    const { jacobian } = constraintData;
    console.log(`Constraint Jacobian rows: ${jacobian?.length || 0}, columns: ${jacobian?.[0]?.length || 0}`);
    
    if (!jacobian || jacobian.length === 0) {
        console.log("No constraints to project against, returning original gradient");
        console.log("======================================================");
        return gradient;
    }
    
    try {
        // Build the inner product matrix
        console.log("Building inner product matrix A_bar");
        const verticesCopy = vertices.map(v => [...v]);
        const result = build_A_bar(
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
        console.log(`Ag dimensions: ${Ag.length}`);
        
        // Solve the saddle point system
        console.log("Solving saddle point system for gradient projection");
        const { x: projectedGradient } = solveSaddlePointSystem(
            A_bar,
            jacobian,
            Ag,
            new Array(jacobian.length).fill(0) // d = 0 for gradient projection
        );
        
        console.log(`Projected gradient length: ${projectedGradient.length} (expected: ${vertices.length * dimension})`);
        if (projectedGradient.length % dimension !== 0) {
            console.error(`CRITICAL: Projected gradient length (${projectedGradient.length}) not divisible by dimension (${dimension})!`);
        }
        
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
 * Works for both 2D and 3D
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
    const dimension = get(config).dimension;
    
    if (!values || values.length === 0) {
        console.log("No constraint values to project, returning original vertices");
        console.log("====================================");
        return vertices.map(v => [...v]);
    }
    
    const initialViolation = Array.isArray(values) ? 
        math.norm(values.filter(v => typeof v === 'number')) : 0;
    console.log(`Initial constraint violation: ${initialViolation}`);
    
    // Working copy of vertices
    let projectedVertices = vertices.map(v => [...v]);
    
    // Settings from config
    const maxIterations = get(config).maxConstraintIterations || 10;
    
    // Track best solution
    let bestVertices = projectedVertices;
    let bestViolation = initialViolation;
    
    console.log(`Will perform at most ${maxIterations} projection iterations`);
    
    // Newton-like iterations as described in Section 5.3.2
    for (let iter = 0; iter < maxIterations; iter++) {
        console.log(`==== Projection iteration ${iter+1}/${maxIterations} ====`);
        
        // Re-evaluate constraint values at current position
        const constraintValues = evaluate(projectedVertices, constraints);
        const constraintNorm = Array.isArray(constraintValues) ? 
            math.norm(constraintValues.filter(v => typeof v === 'number')) : 0;
        
        // Update best solution if this one is better
        if (constraintNorm < bestViolation) {
            bestViolation = constraintNorm;
            bestVertices = projectedVertices.map(v => [...v]);
        }
        
        // Update Jacobian at current position
        console.log("Updating constraint Jacobian at current position");
        const currentJacobian = buildJacobian(projectedVertices, constraints);
        
        try {
            // Build the inner product matrix for the current position
            console.log("Building inner product matrix for constraint projection");
            const result = build_A_bar(
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
            
            // Ensure all constraint values are numeric for the RHS
            const safeConstraintValues = Array.isArray(constraintValues) ? 
                constraintValues.map(v => typeof v === 'number' ? -v : 0) : [];
            
            // Solve the saddle point system
            console.log("Solving saddle point system for constraint projection");
            const { x: correction } = solveSaddlePointSystem(
                A_bar,
                currentJacobian,
                new Array(dimension * projectedVertices.length).fill(0), // b = 0
                safeConstraintValues // d = -Φ(γ̃)
            );
            
            // Apply the correction to vertices based on dimension
            for (let i = 0; i < projectedVertices.length; i++) {
                for (let d = 0; d < dimension; d++) {
                    const correctionVal = correction[i * dimension + d];
                    if (typeof correctionVal === 'number' && isFinite(correctionVal)) {
                        projectedVertices[i][d] += correctionVal;
                    }
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
    const finalViolation = Array.isArray(finalConstraintValues) ? 
        math.norm(finalConstraintValues.filter(v => typeof v === 'number')) : 0;
    console.log(`Final constraint violation: ${finalViolation}`);
    
    // If best solution is significantly better, use that instead
    if (finalViolation > bestViolation * 1.1) {
        console.log(`Final solution is worse than best found, using best found`);
        projectedVertices = bestVertices;
    }
    
    console.log("====================================");
    return projectedVertices;
}