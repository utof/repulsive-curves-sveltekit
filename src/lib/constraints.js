// src/lib/constraints.js
import { get } from 'svelte/store';
import { config } from './stores';
import * as math from 'mathjs';
import { calculateEdgeProperties } from './energyCalculations';

/**
 * Calculates the barycenter of the curve with weighted edges
 *
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {boolean} weighted - Whether to weight by edge lengths
 * @returns {Array} - Barycenter position [x, y]
 */
export function calculateBarycenter(vertices, edges, weighted = true) {
    try {
        console.log("Calculating barycenter");
        // Calculate edge properties if using weighted barycenter
        const { edgeLengths, edgeMidpoints } = weighted ? 
            calculateEdgeProperties(vertices, edges) : { edgeLengths: [], edgeMidpoints: [] };
        
        let barycenter = [0, 0];
        let totalWeight = 0;
        
        if (weighted && edgeLengths && edgeMidpoints) {
            // Weighted by edge lengths as in the paper:
            // Φ_barycenter(γ) := ∑_{I∈E} ℓ_I (x_I - x₀)
            for (let i = 0; i < edges.length; i++) {
                const weight = edgeLengths[i];
                if (!isFinite(weight) || weight <= 0) continue;
                
                barycenter[0] += edgeMidpoints[i][0] * weight;
                barycenter[1] += edgeMidpoints[i][1] * weight;
                totalWeight += weight;
            }
        } else {
            // Simple vertex average (fallback)
            for (const v of vertices) {
                if (!v || v.some(coord => !isFinite(coord))) continue;
                barycenter[0] += v[0];
                barycenter[1] += v[1];
                totalWeight += 1;
            }
        }
        
        if (totalWeight <= 0) {
            console.error("Total weight is zero or negative in barycenter calculation");
            return [0, 0]; // Fallback
        }
        
        barycenter[0] /= totalWeight;
        barycenter[1] /= totalWeight;
        
        return barycenter;
    } catch (error) {
        console.error("Error in calculateBarycenter:", error);
        return [0, 0]; // Fallback
    }
}

/**
 * Evaluates all constraints for the given vertices
 * Returns a flattened array of constraint values
 *
 * @param {Array} vertices - Vertex positions
 * @param {Object} constraints - Constraint configuration
 * @param {Array} edges - Edge connections (optional)
 * @returns {Array} - Constraint values as a flat array
 */
export function evaluateConstraints(vertices, constraints, edges = []) {
    const constraintValues = [];
    
    // Barycenter constraint
    if (constraints.barycenter) {
        console.log("Evaluating barycenter constraint");
        const targetBarycenter = constraints.barycenterTarget || [0, 0];
        const currentBarycenter = calculateBarycenter(vertices, edges, true);
        
        // Constraint is defined as: ∑_{I∈E} ℓ_I (x_I - x₀) = 0
        // This is equivalent to: currentBarycenter - targetBarycenter = 0
        constraintValues.push(currentBarycenter[0] - targetBarycenter[0]);
        constraintValues.push(currentBarycenter[1] - targetBarycenter[1]);
        
        console.log(`Barycenter constraint: current=[${currentBarycenter}], target=[${targetBarycenter}], violation=[${constraintValues}]`);
    }
    
    // Length constraint - currently not implemented but would be added here
    
    return constraintValues;
}

/**
 * Builds the Jacobian matrix for all constraints
 * Each constraint contributes one or more rows to the matrix
 *
 * @param {Array} vertices - Vertex positions
 * @param {Object} constraints - Constraint configuration 
 * @param {Array} edges - Edge connections (optional)
 * @returns {Array} - Jacobian matrix C where C[i][j] = ∂Φᵢ/∂γⱼ
 */
export function buildConstraintJacobian(vertices, constraints, edges = []) {
    const numVertices = vertices.length;
    const jacobian = [];
    
    // Barycenter constraint
    if (constraints.barycenter) {
        console.log("Building Jacobian for barycenter constraint");
        
        // The barycenter constraint has 2 components (x and y)
        // For each component, we need a row in the Jacobian
        
        // For the x-component of barycenter constraint
        const xRow = new Array(numVertices * 2).fill(0);
        
        // For the y-component of barycenter constraint
        const yRow = new Array(numVertices * 2).fill(0);
        
        // Calculate total weight for normalization
        const { edgeLengths } = calculateEdgeProperties(vertices, edges);
        const totalWeight = edgeLengths.reduce((sum, len) => sum + (isFinite(len) ? len : 0), 0);
        
        if (totalWeight <= 0) {
            console.error("Total weight is zero or negative in Jacobian calculation");
            // Add zero rows as fallback
            jacobian.push(xRow);
            jacobian.push(yRow);
            return jacobian;
        }
        
        // We need to determine how the barycenter changes when each vertex moves
        // For a weighted barycenter, each vertex contributes based on its edges
        for (let i = 0; i < edges.length; i++) {
            const [v1, v2] = edges[i];
            const weight = edgeLengths[i] / totalWeight / 2; // Divided by 2 because each vertex contributes half to the midpoint
            
            if (!isFinite(weight)) continue;
            
            // Vertex v1 affects x-coordinate of barycenter
            xRow[v1 * 2] += weight;
            
            // Vertex v1 affects y-coordinate of barycenter
            yRow[v1 * 2 + 1] += weight;
            
            // Vertex v2 affects x-coordinate of barycenter
            xRow[v2 * 2] += weight;
            
            // Vertex v2 affects y-coordinate of barycenter
            yRow[v2 * 2 + 1] += weight;
        }
        
        // Add the rows to the Jacobian
        jacobian.push(xRow);
        jacobian.push(yRow);
        
        console.log(`Barycenter Jacobian: x-row has ${xRow.filter(v => v !== 0).length} non-zero entries`);
        console.log(`Barycenter Jacobian: y-row has ${yRow.filter(v => v !== 0).length} non-zero entries`);
    }
    
    // Length constraint Jacobian would be added here
    
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
    try {
        console.log("Creating constraint data");
        
        // Validate edges parameter
        if (!edges || !Array.isArray(edges)) {
            console.error("Invalid edges parameter:", edges);
            edges = [];
        }
        
        // Evaluate constraints
        const values = evaluateConstraints(vertices, constraints, edges);
        
        // Build Jacobian
        const jacobian = buildConstraintJacobian(vertices, constraints, edges);
        
        return {
            values,
            jacobian,
            evaluate: (v, c) => evaluateConstraints(v, c, edges),
            buildJacobian: (v, c) => buildConstraintJacobian(v, c, edges)
        };
    } catch (error) {
        console.error("Error creating constraint data:", error);
        return {
            values: [],
            jacobian: [],
            evaluate: () => [],
            buildJacobian: () => []
        };
    }
}


// Export relevant constraint types for reference
export const ConstraintTypes = {
    BARYCENTER: 'barycenter',
    LENGTH: 'length',
    EDGE_LENGTH: 'edge_length',
    POINT: 'point',
    SURFACE: 'surface',
    TANGENT: 'tangent'
};