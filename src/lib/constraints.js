// src/lib/constraints.js
import { get } from 'svelte/store';
import { config, initialTotalLength } from './stores';
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
 * Calculate the total length of a curve
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @returns {number} Total curve length
 */
export function calculateTotalLength(vertices, edges) {
    try {
        let totalLength = 0;
        for (const [v1, v2] of edges) {
            const dx = vertices[v2][0] - vertices[v1][0];
            const dy = vertices[v2][1] - vertices[v1][1];
            totalLength += Math.sqrt(dx*dx + dy*dy);
        }
        return totalLength;
    } catch (error) {
        console.error("Error calculating total length:", error);
        return 0;
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
        
        // Get scaling factor for better numerical conditioning
        const scalingFactor = get(config).barycenterScaling || 0.01;
        
        // Constraint is defined as: ∑_{I∈E} ℓ_I (x_I - x₀) = 0
        // This is equivalent to: (currentBarycenter - targetBarycenter) * scalingFactor = 0
        constraintValues.push((currentBarycenter[0] - targetBarycenter[0]) * scalingFactor);
        constraintValues.push((currentBarycenter[1] - targetBarycenter[1]) * scalingFactor);
        
        console.log(`Barycenter constraint: current=[${currentBarycenter}], target=[${targetBarycenter}], scaled violation=[${constraintValues.slice(-2)}]`);
    }
    
    // Length constraint
    if (constraints.length) {
        console.log("Evaluating length constraint");
        
        // Target length can be a direct value or a percentage of initial length
        let targetLength = constraints.lengthTarget || 0;
        
        // Check if target should be calculated from percentage
        if (constraints.lengthPercentage) {
            const initialLength = get(initialTotalLength);
            targetLength = initialLength * (constraints.lengthPercentage / 100);
            console.log(`Using ${constraints.lengthPercentage}% of initial length (${initialLength.toFixed(2)}) as target: ${targetLength.toFixed(2)}`);
        }
        
        const currentLength = calculateTotalLength(vertices, edges);
        
        // Get scaling factor for better numerical conditioning
        const scalingFactor = get(config).lengthScaling || 0.001;
        console.log(`Using length stabilization factor: ${scalingFactor}`);
        
        // Constraint is defined as: (L^0 - ∑_{I∈E} ℓ_I) * scalingFactor = 0
        constraintValues.push((targetLength - currentLength) * scalingFactor);
        
        console.log(`Length constraint: current=${currentLength.toFixed(2)}, target=${targetLength.toFixed(2)}, scaled violation=${((targetLength - currentLength) * scalingFactor).toFixed(4)}`);
    }
    
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
        
        // Calculate edge properties for weighted contribution calculation
        const { edgeLengths } = calculateEdgeProperties(vertices, edges);
        const totalWeight = edgeLengths.reduce((sum, len) => sum + (isFinite(len) ? len : 0), 0);
        
        if (totalWeight <= 0) {
            console.error("Total weight is zero or negative in Jacobian calculation");
            // Add zero rows as fallback
            jacobian.push(xRow);
            jacobian.push(yRow);
            return jacobian;
        }
        
        // Scale factor for the Jacobian to improve numerical stability
        const stabilizationFactor = get(config).barycenterScaling || 0.01;
        
        console.log(`Using barycenter stabilization factor: ${stabilizationFactor}`);
        
        // For each edge, compute how vertices affect the weighted barycenter
        for (let i = 0; i < edges.length; i++) {
            const [v1, v2] = edges[i];
            // Weight is normalized by total weight for better conditioning
            // Each vertex contributes half to the midpoint
            const weight = edgeLengths[i] / totalWeight / 2 * stabilizationFactor;
            
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
        
        // Log max values for diagnostics
        const xRowMax = Math.max(...xRow.map(Math.abs));
        const yRowMax = Math.max(...yRow.map(Math.abs));
        console.log(`Barycenter Jacobian: x-row max value=${xRowMax}, y-row max value=${yRowMax}`);
    }
    
    // Length constraint Jacobian
    if (constraints.length) {
        console.log("Building Jacobian for length constraint");
        
        // One row for the length constraint
        const lengthRow = new Array(numVertices * 2).fill(0);
        
        // Scale factor for the Jacobian to improve numerical stability
        const stabilizationFactor = get(config).lengthScaling || 0.001;
        
        console.log(`Using length stabilization factor: ${stabilizationFactor}`);
        
        // For each edge, calculate derivatives of length with respect to vertices
        for (const [v1, v2] of edges) {
            const dx = vertices[v2][0] - vertices[v1][0];
            const dy = vertices[v2][1] - vertices[v1][1];
            const length = Math.sqrt(dx*dx + dy*dy);
            
            if (!isFinite(length) || length <= 0) continue;
            
            // Derivatives of length with respect to each vertex component
            // d(length)/d(v1.x) = -dx/length
            // d(length)/d(v1.y) = -dy/length
            // d(length)/d(v2.x) = dx/length
            // d(length)/d(v2.y) = dy/length
            
            // Apply stabilization factor to improve numerical conditioning
            const factor = stabilizationFactor / length;
            
            // v1 derivatives (negative because moving v1 away from v2 increases length)
            lengthRow[v1 * 2] -= dx * factor;
            lengthRow[v1 * 2 + 1] -= dy * factor;
            
            // v2 derivatives
            lengthRow[v2 * 2] += dx * factor;
            lengthRow[v2 * 2 + 1] += dy * factor;
        }
        
        // Add row to the Jacobian
        jacobian.push(lengthRow);
        
        // Log max value for diagnostics
        const lengthRowMax = Math.max(...lengthRow.map(Math.abs));
        console.log(`Length Jacobian: max value=${lengthRowMax}, non-zero entries=${lengthRow.filter(v => v !== 0).length}`);
    }
    
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
        
        // Return object with constraint data and functions to re-evaluate
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