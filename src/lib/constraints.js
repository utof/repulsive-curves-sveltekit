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
 * @returns {Array} - Barycenter position [x, y]
 */
export function calculateBarycenter(vertices, edges, weighted = true) {
    console.log("==== CALCULATING BARYCENTER ====");
    // Calculate edge properties if using weighted barycenter
    const edgeLengths = [];
    const edgeMidpoints = [];
    
    if (weighted) {
        for (const [v1, v2] of edges) {
            const dx = vertices[v2][0] - vertices[v1][0];
            const dy = vertices[v2][1] - vertices[v1][1];
            const length = Math.sqrt(dx*dx + dy*dy);
            
            edgeLengths.push(length);
            edgeMidpoints.push([
                (vertices[v1][0] + vertices[v2][0]) / 2,
                (vertices[v1][1] + vertices[v2][1]) / 2
            ]);
        }
    }
    
    let barycenter = [0, 0];
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
            
            barycenter[0] += edgeMidpoints[i][0] * weight;
            barycenter[1] += edgeMidpoints[i][1] * weight;
            totalWeight += weight;
        }
    } else {
        // Simple vertex average (fallback)
        console.log("Using simple vertex average (not weighted)");
        for (const v of vertices) {
            if (!v || v.some(coord => !isFinite(coord))) continue;
            barycenter[0] += v[0];
            barycenter[1] += v[1];
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
    barycenter[0] /= totalWeight;
    barycenter[1] /= totalWeight;
    
    console.log(`Calculated barycenter: [${barycenter[0].toFixed(4)}, ${barycenter[1].toFixed(4)}]`);
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
    
    for (let i = 0; i < edges.length; i++) {
        const [v1, v2] = edges[i];
        const dx = vertices[v2][0] - vertices[v1][0];
        const dy = vertices[v2][1] - vertices[v1][1];
        const length = Math.sqrt(dx*dx + dy*dy);
        
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
 * @returns {Array} - [x-component, y-component] of constraint values
 */
export function evaluateBarycenterConstraint(vertices, edges, targetBarycenter, scalingFactor) {
    console.log("==== EVALUATING BARYCENTER CONSTRAINT ====");
    console.log(`Target barycenter: [${targetBarycenter[0].toFixed(4)}, ${targetBarycenter[1].toFixed(4)}]`);
    
    // Get current barycenter
    const currentBarycenter = calculateBarycenter(vertices, edges, true);
    
    // Get scaling factor for better numerical conditioning
    if (scalingFactor === undefined) {
        scalingFactor = ConstraintNumerics.USE_STABILITY_SCALING ? 
            get(config).barycenterScaling || ConstraintNumerics.DEFAULT_BARYCENTER_SCALING : 1.0;
    }
    
    console.log(`Using barycenter scaling factor: ${scalingFactor}`);
    
    // Constraint value: (current - target) * scaling
    const constraintX = (currentBarycenter[0] - targetBarycenter[0]) * scalingFactor;
    const constraintY = (currentBarycenter[1] - targetBarycenter[1]) * scalingFactor;
    
    console.log(`Barycenter constraint values: [${constraintX.toFixed(6)}, ${constraintY.toFixed(6)}]`);
    console.log("=========================================");
    
    return [constraintX, constraintY];
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
        const dx = vertices[v2][0] - vertices[v1][0];
        const dy = vertices[v2][1] - vertices[v1][1];
        const currentLength = Math.sqrt(dx*dx + dy*dy);
        
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
    
    // Barycenter constraint
    if (constraints.barycenter) {
        console.log("Processing barycenter constraint");
        const targetBarycenter = constraints.barycenterTarget || [0, 0];
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