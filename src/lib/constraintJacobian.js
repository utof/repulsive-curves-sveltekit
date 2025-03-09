// src/lib/constraintJacobian.js
import { get } from 'svelte/store';
import { config } from './stores';
import { 
    calculateBarycenter,
    calculateTotalLength,
    evaluateConstraints,
    ConstraintNumerics
} from './constraints';

/**
 * Builds the constraint Jacobian matrix using finite differences
 * More general than analytical formulas and works for any constraint
 * Works for both 2D and 3D
 * 
 * @param {Array} vertices - Vertex positions
 * @param {Object} constraints - Constraint configuration
 * @param {Array} edges - Edge connections
 * @returns {Array} - Jacobian matrix C where C[i][j] = ∂Φᵢ/∂γⱼ
 */
export function buildConstraintJacobian(vertices, constraints, edges) {
    console.log("==== BUILDING CONSTRAINT JACOBIAN USING FINITE DIFFERENCES ====");
    
    // Get current dimension
    const dimension = get(config).dimension;
    const numVertices = vertices.length;
    
    // Step size for finite differences
    const h = get(config).finiteDiffH || 1e-7;
    console.log(`Using finite difference step size h = ${h}`);
    
    // First, evaluate the constraints at the current position
    const baseConstraintValues = evaluateConstraints(vertices, constraints, edges);
    console.log(`Base constraint values: [${baseConstraintValues.map(v => typeof v === 'number' ? v.toFixed(6) : v).join(', ')}]`);
    
    const numConstraints = baseConstraintValues.length;
    if (numConstraints === 0) {
        console.log("No constraints to evaluate");
        return [];
    }
    
    // Initialize the Jacobian matrix with zeros
    // Each row corresponds to one constraint
    // Each column corresponds to one vertex coordinate
    const jacobian = [];
    for (let i = 0; i < numConstraints; i++) {
        jacobian.push(new Array(numVertices * dimension).fill(0));
    }
    
    // For each vertex and each dimension, compute partial derivatives
    for (let vertexIdx = 0; vertexIdx < numVertices; vertexIdx++) {
        console.log(`Computing derivatives for vertex ${vertexIdx}`);
        
        for (let dim = 0; dim < dimension; dim++) {
            // Create a copy of the vertices
            const perturbedVertices = vertices.map(v => [...v]);
            
            // Perturb the current vertex coordinate
            const originalValue = perturbedVertices[vertexIdx][dim];
            perturbedVertices[vertexIdx][dim] += h;
            
            // Evaluate the constraints with the perturbed position
            const perturbedConstraintValues = evaluateConstraints(
                perturbedVertices, 
                constraints, 
                edges
            );
            
            // Calculate the partial derivatives using finite differences
            for (let constraintIdx = 0; constraintIdx < numConstraints; constraintIdx++) {
                const baseValue = baseConstraintValues[constraintIdx];
                const perturbedValue = perturbedConstraintValues[constraintIdx];
                
                // Skip if values are not numeric
                if (typeof baseValue !== 'number' || typeof perturbedValue !== 'number') {
                    console.warn(`Non-numeric constraint values detected at constraint ${constraintIdx}`);
                    continue;
                }
                
                // Calculate derivative: (f(x+h) - f(x)) / h
                const derivative = (perturbedValue - baseValue) / h;
                
                // Store in Jacobian
                const columnIndex = vertexIdx * dimension + dim;
                jacobian[constraintIdx][columnIndex] = derivative;
            }
            
            // Restore the original value
            perturbedVertices[vertexIdx][dim] = originalValue;
        }
    }
    
    // Log statistics about the Jacobian
    let nonZeroCount = 0;
    let maxValue = 0;
    let minAbsNonZero = Infinity;
    
    for (const row of jacobian) {
        for (const value of row) {
            if (value !== 0) {
                nonZeroCount++;
                maxValue = Math.max(maxValue, Math.abs(value));
                minAbsNonZero = Math.min(minAbsNonZero, Math.abs(value));
            }
        }
    }
    
    console.log(`Jacobian stats: ${nonZeroCount} non-zero entries, max abs value: ${maxValue.toExponential(4)}, min abs non-zero: ${minAbsNonZero.toExponential(4)}`);
    console.log("===========================================================");
    
    return jacobian;
}

/**
 * Creates a constraint data object containing all information needed for constraint handling
 *
 * @param {Array} vertices - Vertex positions
 * @param {Object} constraints - Constraint configuration
 * @param {Array} edges - Edge connections
 * @returns {Object} - Constraint data including values and Jacobian
 */
export function createConstraintData(vertices, constraints, edges) {
    console.log("==== CREATING CONSTRAINT DATA ====");
    
    // Evaluate constraints
    const values = evaluateConstraints(vertices, constraints, edges);
    
    // Build Jacobian
    const jacobian = buildConstraintJacobian(vertices, constraints, edges);
    
    // Return object with constraint data and functions to re-evaluate
    const constraintData = {
        values,
        jacobian,
        evaluate: (v, c) => evaluateConstraints(v, c, edges),
        buildJacobian: (v, c) => buildConstraintJacobian(v, c, edges)
    };
    
    console.log(`Created constraint data with ${values.length} values and ${jacobian.length} Jacobian rows`);
    console.log("================================");
    
    return constraintData;
}