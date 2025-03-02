// src/lib/constraints.js
import { get } from 'svelte/store';
import { config } from './stores';
import * as math from 'mathjs';
import { calculateEdgeProperties } from './energyCalculations';

/**
 * Enforces barycenter constraint on vertices
 * 
 * According to the paper, the barycenter constraint is defined as:
 * Φ_barycenter(γ) := ∑_{I∈E} ℓ_I (x_I - x₀) = 0
 * 
 * Where:
 * - γ is the curve (vertices)
 * - E is the set of edges
 * - ℓ_I is the length of edge I
 * - x_I is the midpoint of edge I
 * - x₀ is the target barycenter position
 * 
 * This constraint ensures that the weighted average of edge midpoints 
 * (weighted by edge length) stays at the specified target position x₀.
 *
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Object} options - Barycenter options
 * @returns {Array} - Vertices with enforced barycenter
 */
export function enforceBarycenter(vertices, edges, options = {}) {
    try {
        console.log("Enforcing barycenter constraint");
        
        const {
            targetBarycenter = [0, 0],
            weighted = true
        } = options;
        
        // Make a copy of vertices to avoid mutating the original
        const projectedVertices = vertices.map(v => [...v]);
        
        // Calculate edge properties - lengths and midpoints
        const { edgeLengths, edgeMidpoints } = calculateEdgeProperties(projectedVertices, edges);
        
        if (!edgeLengths || !edgeMidpoints) {
            console.error("Failed to calculate edge properties", { edgeLengths, edgeMidpoints });
            return projectedVertices;
        }
        
        // Current barycenter calculation
        let currentBarycenter = [0, 0];
        let totalWeight = 0;
        
        if (weighted) {
            // Use edge-length-weighted barycenter as defined in the paper:
            // Barycenter = ∑_{I∈E} ℓ_I × x_I / ∑_{I∈E} ℓ_I
            for (let i = 0; i < edges.length; i++) {
                const weight = edgeLengths[i];
                
                if (weight <= 0 || !isFinite(weight)) {
                    console.warn(`Invalid edge length for edge ${i}: ${weight}`);
                    continue;
                }
                
                if (!edgeMidpoints[i] || edgeMidpoints[i].some(v => !isFinite(v))) {
                    console.warn(`Invalid midpoint for edge ${i}: ${edgeMidpoints[i]}`);
                    continue;
                }
                
                currentBarycenter[0] += edgeMidpoints[i][0] * weight;
                currentBarycenter[1] += edgeMidpoints[i][1] * weight;
                totalWeight += weight;
            }
        } else {
            // Alternative: Use simple vertex average (not from the paper)
            for (const v of projectedVertices) {
                if (!v || v.some(coord => !isFinite(coord))) {
                    console.warn("Skipping invalid vertex:", v);
                    continue;
                }
                currentBarycenter[0] += v[0];
                currentBarycenter[1] += v[1];
            }
            totalWeight = projectedVertices.length;
        }
        
        if (totalWeight <= 0) {
            console.error("Total weight is zero or negative, cannot enforce barycenter", { totalWeight });
            return projectedVertices;
        }
        
        currentBarycenter[0] /= totalWeight;
        currentBarycenter[1] /= totalWeight;
        
        // Calculate translation needed
        const dx = targetBarycenter[0] - currentBarycenter[0];
        const dy = targetBarycenter[1] - currentBarycenter[1];
        
        if (!isFinite(dx) || !isFinite(dy)) {
            console.error("Invalid translation vector", { dx, dy });
            return projectedVertices;
        }
        
        console.log(`Barycenter translation: (${dx.toFixed(4)}, ${dy.toFixed(4)})`);
        
        // Apply translation to all vertices
        for (let i = 0; i < projectedVertices.length; i++) {
            projectedVertices[i][0] += dx;
            projectedVertices[i][1] += dy;
        }
        
        return projectedVertices;
    } catch (error) {
        console.error("Error in enforceBarycenter:", error);
        // Return the original vertices if there's an error
        return vertices.map(v => [...v]);
    }
}

/**
 * Applies constraints in sequence
 * Currently only the barycenter constraint is implemented
 * 
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Object} constraints - Constraint configuration
 * @param {Object} additionalData - Additional data for constraints
 * @returns {Array} - Vertices after applying all constraints
 */
export function applyConstraints(vertices, edges, constraints, additionalData = {}) {
    let result = vertices.map(v => [...v]);
    
    // Apply barycenter constraint if enabled
    if (constraints.barycenter) {
        console.log(`Applying barycenter constraint: ${JSON.stringify(constraints.barycenterTarget)}`);
        result = enforceBarycenter(result, edges, {
            targetBarycenter: constraints.barycenterTarget || [0, 0],
            weighted: true
        });
    }
    
    return result;
}

/**
 * Projects gradient onto the tangent space of the barycenter constraint
 * 
 * Section 5.3.1 in the paper describes this as finding a descent direction
 * that is tangent to the constraint set.
 * 
 * For the barycenter constraint, this means removing any component of the 
 * gradient that would change the barycenter.
 *
 * @param {Array} gradient - Gradient to project
 * @param {Array} vertices - Current vertex positions
 * @param {Array} edges - Edge connections
 * @param {Object} constraints - Active constraints configuration
 * @returns {Array} - Projected gradient
 */
export function projectGradientOntoConstraints(gradient, vertices, edges, constraints) {
    try {
        // If no barycenter constraint, return the original gradient
        if (!constraints.barycenter) {
            return gradient.map(g => [...g]);
        }
        
        console.log("Projecting gradient onto barycenter constraint tangent space");
        
        // Make a copy of the gradient
        let projectedGradient = gradient.map(g => [...g]);
        
        // For barycenter constraint, we need to remove the translation component
        // Calculate average gradient (translation component)
        const avgGradient = [0, 0];
        let validGradientCount = 0;
        
        for (const g of projectedGradient) {
            if (!g || g.some(v => !isFinite(v))) {
                console.warn("Skipping invalid gradient component:", g);
                continue;
            }
            avgGradient[0] += g[0];
            avgGradient[1] += g[1];
            validGradientCount++;
        }
        
        if (validGradientCount === 0) {
            console.error("No valid gradient components found");
            return gradient.map(g => [...g]); // Return original gradient
        }
        
        avgGradient[0] /= validGradientCount;
        avgGradient[1] /= validGradientCount;
        
        console.log(`Translation component to remove: (${avgGradient[0].toFixed(6)}, ${avgGradient[1].toFixed(6)})`);
        
        // Subtract translation component from all gradients
        for (let i = 0; i < projectedGradient.length; i++) {
            projectedGradient[i][0] -= avgGradient[0];
            projectedGradient[i][1] -= avgGradient[1];
        }
        
        return projectedGradient;
    } catch (error) {
        console.error("Error in projectGradientOntoConstraints:", error);
        return gradient.map(g => [...g]); // Return original gradient on error
    }
}