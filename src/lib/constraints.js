// src/lib/constraints.js
import { get } from 'svelte/store';
import { config, initialTotalLength } from './stores';

// Numerical stability toggles and parameters
export const ConstraintNumerics = {
    // Toggle numerical stability measures
    USE_STABILITY_SCALING: false,     // Use scaling factors for better conditioning
    USE_SAFE_DENOMINATORS: false,     // Avoid division by zero
    
    // Parameter values (only used if stability measures are enabled)
    MIN_LENGTH: 1e-10,              // Minimum edge length to consider
    MIN_TOTAL_WEIGHT: 1e-8,         // Minimum total weight in barycenter
    DEFAULT_BARYCENTER_SCALING: 1.0, // Scaling factor for barycenter constraint
    DEFAULT_LENGTH_SCALING: 1.0,     // Scaling factor for length constraint
};

/**
 * Calculates the barycenter of the curve with weighted edges
 * Following the formulation in the paper: 
 * Φ_barycenter(y) := ∑_{I∈E} ℓ_I (x_I - x₀)
 *
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {boolean} weighted - Whether to weight by edge lengths
 * @returns {Array} - Barycenter position [x, y] or [x, y, z] for 3D
 */
export function calculateBarycenter(vertices, edges, weighted = true) {
    console.log("==== CALCULATING BARYCENTER ====");
    const dimension = get(config).dimension;
    
    // Calculate edge properties if using weighted barycenter
    const edgeLengths = [];
    const edgeMidpoints = [];
    
    if (weighted) {
        for (const [v1, v2] of edges) {
            const v1Pos = vertices[v1];
            const v2Pos = vertices[v2];
            
            // // Safely check vertices
            // if (!v1Pos || !v2Pos) continue;
            
            // Calculate edge vector and length
            const diff = [];
            for (let i = 0; i < dimension; i++) {
                diff[i] = v2Pos[i] - v1Pos[i];
            }
            const length = Math.sqrt(diff.reduce((sum, val) => sum + val * val, 0));
            
            // Calculate midpoint
            const midpoint = [];
            for (let i = 0; i < dimension; i++) {
                midpoint[i] = (v1Pos[i] + v2Pos[i]) / 2;
            }
            
            edgeLengths.push(length);
            edgeMidpoints.push(midpoint);
        }
    }
    
    let barycenter = new Array(dimension).fill(0);
    let totalWeight = 0;
    
    if (weighted && edgeLengths.length > 0) {
        // Weighted by edge lengths as in the paper
        console.log("Using weighted barycenter calculation");
        for (let i = 0; i < edges.length; i++) {
            const weight = edgeLengths[i];
            
            // Skip invalid lengths
            if (!isFinite(weight) || 
                (ConstraintNumerics.USE_SAFE_DENOMINATORS && weight <= ConstraintNumerics.MIN_LENGTH)) {
                console.log(`Skipping edge ${i} with invalid length ${weight}`);
                continue;
            }
            
            // Add weighted contribution for each dimension
            for (let d = 0; d < dimension; d++) {
                barycenter[d] += edgeMidpoints[i][d] * weight;
            }
            totalWeight += weight;
        }
    } else {
        // Simple vertex average (fallback)
        console.log("Using simple vertex average (not weighted)");
        for (const v of vertices) {
            // if (!v || v.some(coord => !isFinite(coord))) continue;
            
            for (let d = 0; d < dimension; d++) {
                barycenter[d] += v[d];
            }
            totalWeight += 1;
        }
    }
    
    // Check total weight
    if (totalWeight <= 0 || !isFinite(totalWeight)) {
        console.warn(`Invalid total weight in barycenter calculation: ${totalWeight}`);
        if (ConstraintNumerics.USE_SAFE_DENOMINATORS) {
            totalWeight = ConstraintNumerics.MIN_TOTAL_WEIGHT;
            console.log(`Using minimum weight value: ${totalWeight}`);
        }
    }
    
    // Compute final barycenter
    for (let d = 0; d < dimension; d++) {
        barycenter[d] /= totalWeight;
    }
    
    console.log(`Calculated barycenter: [${barycenter.map(v => v.toFixed(4)).join(', ')}]`);
    console.log("================================");
    
    return barycenter;
}

/**
 * Calculate the total length of a curve
 * L = ∑_{I∈E} ℓ_I
 * 
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @returns {number} Total curve length
 */
export function calculateTotalLength(vertices, edges) {
    console.log("==== CALCULATING TOTAL LENGTH ====");
    let totalLength = 0;
    const dimension = get(config).dimension;
    
    for (let i = 0; i < edges.length; i++) {
        const [v1, v2] = edges[i];
        const v1Pos = vertices[v1];
        const v2Pos = vertices[v2];
        
        // if (!v1Pos || !v2Pos) continue;
        
        // Calculate squared distance
        let squaredDist = 0;
        for (let d = 0; d < dimension; d++) {
            const diff = v2Pos[d] - v1Pos[d];
            squaredDist += diff * diff;
        }
        const length = Math.sqrt(squaredDist);
        
        if (isFinite(length)) {
            totalLength += length;
            console.log(`Edge ${i}: [${v1}, ${v2}], length = ${length.toFixed(4)}`);
        } else {
            console.warn(`Edge ${i}: [${v1}, ${v2}] has invalid length`);
        }
    }
    
    console.log(`Total length: ${totalLength.toFixed(4)}`);
    console.log("================================");
    
    return totalLength;
}

/**
 * Evaluates barycenter constraint 
 * Constraint is: Φ_barycenter(γ) := ∑_{I∈E} ℓ_I (x_I - x₀) = 0
 * 
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Array} targetBarycenter - Target barycenter position
 * @param {number} scalingFactor - Optional scaling for numerical stability
 * @returns {Array} - [x, y, z?] Constraint values, one per dimension
 */
export function evaluateBarycenterConstraint(vertices, edges, targetBarycenter, scalingFactor) {
    console.log("==== EVALUATING BARYCENTER CONSTRAINT ====");
    console.log(`Target barycenter: [${targetBarycenter.map(v => v.toFixed(4)).join(', ')}]`);
    
    const dimension = get(config).dimension;
    
    // Get current barycenter
    const currentBarycenter = calculateBarycenter(vertices, edges, true);
    
    // Get scaling factor for better numerical conditioning
    if (scalingFactor === undefined) {
        scalingFactor = ConstraintNumerics.USE_STABILITY_SCALING ? 
            get(config).barycenterScaling || ConstraintNumerics.DEFAULT_BARYCENTER_SCALING : 1.0;
    }
    
    console.log(`Using barycenter scaling factor: ${scalingFactor}`);
    
    // Constraint values: (current - target) * scaling
    const constraintValues = [];
    for (let d = 0; d < dimension; d++) {
        // If target doesn't specify this dimension (e.g., Z for a 2D target in 3D space),
        // use the current value (i.e., don't constrain this dimension)
        const targetValue = d < targetBarycenter.length ? targetBarycenter[d] : currentBarycenter[d];
        const constraintValue = (currentBarycenter[d] - targetValue) * scalingFactor;
        constraintValues.push(constraintValue);
    }
    
    console.log(`Barycenter constraint values: [${constraintValues.map(v => v.toFixed(6)).join(', ')}]`);
    console.log("=========================================");
    
    return constraintValues;
}

/**
 * Evaluates length constraint
 * Constraint is: Φ_length(γ) := L^0 - ∑_{I∈E} ℓ_I = 0
 * 
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Object} lengthConstraint - Length constraint settings
 * @param {number} scalingFactor - Optional scaling for numerical stability
 * @returns {number} - Constraint value
 */
export function evaluateLengthConstraint(vertices, edges, lengthConstraint, scalingFactor) {
    console.log("==== EVALUATING LENGTH CONSTRAINT ====");
    
    // Get current length
    const currentLength = calculateTotalLength(vertices, edges);
    
    // Determine target length
    let targetLength = lengthConstraint.lengthTarget || 0;
    
    // Check if target should be calculated from percentage
    if (lengthConstraint.lengthPercentage) {
        const initialLength = get(initialTotalLength);
        targetLength = initialLength * (lengthConstraint.lengthPercentage / 100);
        console.log(`Using ${lengthConstraint.lengthPercentage}% of initial length (${initialLength.toFixed(2)}) as target: ${targetLength.toFixed(2)}`);
    } else {
        console.log(`Using fixed target length: ${targetLength.toFixed(2)}`);
    }
    
    // Get scaling factor for better numerical conditioning
    if (scalingFactor === undefined) {
        scalingFactor = ConstraintNumerics.USE_STABILITY_SCALING ? 
            get(config).lengthScaling || ConstraintNumerics.DEFAULT_LENGTH_SCALING : 1.0;
    }
    
    console.log(`Using length scaling factor: ${scalingFactor}`);
    
    // Constraint is: (L^0 - ∑_{I∈E} ℓ_I) * scalingFactor = 0
    const constraintValue = (targetLength - currentLength) * scalingFactor;
    
    console.log(`Length constraint: target=${targetLength.toFixed(2)}, current=${currentLength.toFixed(2)}, scaled violation=${constraintValue.toFixed(6)}`);
    console.log("====================================");
    
    return constraintValue;
}

/**
 * Evaluates edge length constraints
 * Constraint is: Φ_edge_length,I(γ) := ℓ_I^0 - ℓ_I = 0
 * Where ℓ_I^0 is the target length for each edge
 * 
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Array} targetLengths - Target length for each edge
 * @param {number} scalingFactor - Optional scaling for numerical stability
 * @returns {Array} - Array of constraint values, one per edge
 */
export function evaluateEdgeLengthConstraints(vertices, edges, targetLengths, scalingFactor) {
    console.log("==== EVALUATING EDGE LENGTH CONSTRAINTS ====");
    const constraintValues = [];
    const dimension = get(config).dimension;
    
    if (!targetLengths || targetLengths.length === 0) {
        console.warn("No target lengths provided for edge length constraints");
        return [];
    }
    
    // Get scaling factor for better numerical conditioning
    if (scalingFactor === undefined) {
        scalingFactor = ConstraintNumerics.USE_STABILITY_SCALING ? 
            get(config).edgeLengthScaling || ConstraintNumerics.DEFAULT_LENGTH_SCALING : 1.0;
    }
    
    console.log(`Using edge length scaling factor: ${scalingFactor}`);
    
    for (let i = 0; i < edges.length; i++) {
        const [v1, v2] = edges[i];
        const v1Pos = vertices[v1];
        const v2Pos = vertices[v2];
        
        // Skip if vertices don't exist
        if (!v1Pos || !v2Pos) {
            constraintValues.push(0);
            console.error(`Edge ${i} has missing vertices`);
            continue;
        }
        
        // Calculate edge length
        let squaredDist = 0;
        for (let d = 0; d < dimension; d++) {
            const diff = v2Pos[d] - v1Pos[d];
            squaredDist += diff * diff;
        }
        const currentLength = Math.sqrt(squaredDist);
        
        // Use the target length for this edge
        const targetLength = targetLengths[i % targetLengths.length];
        
        // Constraint is: (ℓ_I^0 - ℓ_I) * scalingFactor = 0
        const constraintValue = (targetLength - currentLength) * scalingFactor;
        constraintValues.push(constraintValue);
        
        console.log(`Edge ${i}: target=${targetLength.toFixed(4)}, current=${currentLength.toFixed(4)}, violation=${constraintValue.toFixed(6)}`);
    }
    
    console.log("===========================================");
    return constraintValues;
}


/**
 * Evaluates all constraints for the given vertices
 * Returns a flattened array of constraint values
 *
 * @param {Array} vertices - Vertex positions
 * @param {Object} constraints - Constraint configuration
 * @param {Array} edges - Edge connections
 * @returns {Array} - Constraint values as a flat array
 */
export function evaluateConstraints(vertices, constraints, edges) {
    console.log("==== EVALUATING ALL CONSTRAINTS ====");
    const constraintValues = [];
    const dimension = get(config).dimension;
    
    // Barycenter constraint
    if (constraints.barycenter) {
        console.log("Processing barycenter constraint");
        const targetBarycenter = constraints.barycenterTarget || new Array(dimension).fill(0);
        const barycenterValues = evaluateBarycenterConstraint(
            vertices, 
            edges, 
            targetBarycenter
        );
        
        constraintValues.push(...barycenterValues);
    }
    
    // Length constraint
    if (constraints.length) {
        console.log("Processing length constraint");
        const lengthConstraint = {
            lengthTarget: constraints.lengthTarget,
            lengthPercentage: constraints.lengthPercentage
        };
        
        const lengthValue = evaluateLengthConstraint(
            vertices,
            edges,
            lengthConstraint
        );
        
        constraintValues.push(lengthValue);
    }
    
    if (constraints.edgeLength && constraints.edgeLengthTargets) {
        console.log("Processing edge length constraints");
        
        const edgeLengthValues = evaluateEdgeLengthConstraints(
            vertices,
            edges,
            constraints.edgeLengthTargets
        );
        
        constraintValues.push(...edgeLengthValues);
        console.log(`Added ${edgeLengthValues.length} edge length constraint values`);
    }
    console.log(`Total constraint values: [${constraintValues.map(v => v.toFixed(6)).join(', ')}]`);
    console.log("==================================");
    
    return constraintValues;
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