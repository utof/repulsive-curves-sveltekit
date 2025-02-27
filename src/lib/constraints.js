// src/lib/constraints.js
import { get } from 'svelte/store';
import { config } from './stores';
import * as math from 'mathjs';
import { calculateEdgeProperties } from './energyCalculations';

/**
 * Projects vertices to maintain original edge lengths
 * Simple but robust approach that directly scales edges toward target lengths
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Array} initialEdgeLengths - Original edge lengths to maintain
 * @param {Number} maxIterations - Maximum projection iterations
 * @param {Number} tolerance - Error tolerance for projection
 * @returns {Array} - Projected vertex positions
 */
export function projectConstraints(
    vertices,
    edges,
    initialEdgeLengths,
    maxIterations = get(config).maxConstraintIterations || 10,
    tolerance = get(config).constraintTolerance || 1e-2
) {
    const projectedVertices = vertices.map(v => [...v]);
    
    for (let iter = 0; iter < maxIterations; iter++) {
        let maxError = 0;
        
        for (let i = 0; i < edges.length; i++) {
            const [v1Idx, v2Idx] = edges[i];
            const v1 = projectedVertices[v1Idx];
            const v2 = projectedVertices[v2Idx];
            
            const dx = v2[0] - v1[0];
            const dy = v2[1] - v1[1];
            const currentLength = Math.sqrt(dx * dx + dy * dy);
            const targetLength = initialEdgeLengths[i];
            
            if (currentLength < 1e-8) continue; // Skip degenerate edges
            
            const error = Math.abs(currentLength - targetLength) / targetLength;
            maxError = Math.max(maxError, error);
            
            // Scale edge to target length
            if (error > tolerance) {
                const scale = targetLength / currentLength;
                const midX = (v1[0] + v2[0]) / 2;
                const midY = (v1[1] + v2[1]) / 2;
                
                // Scale from midpoint
                projectedVertices[v1Idx][0] = midX - (dx / 2) * scale;
                projectedVertices[v1Idx][1] = midY - (dy / 2) * scale;
                projectedVertices[v2Idx][0] = midX + (dx / 2) * scale;
                projectedVertices[v2Idx][1] = midY + (dy / 2) * scale;
            }
        }
        
        if (maxError < tolerance) {
            console.log(`Length constraint satisfied after ${iter + 1} iterations with error ${maxError}`);
            break;
        }
    }
    
    return projectedVertices;
}

/**
 * Enforces barycenter constraint on vertices
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Object} options - Barycenter options
 * @returns {Array} - Vertices with enforced barycenter
 */
export function enforceBarycenter(vertices, edges, options = {}) {
    const {
        targetBarycenter = [0, 0],
        weighted = true
    } = options;
    
    const projectedVertices = vertices.map(v => [...v]);
    
    // Calculate current barycenter
    let currentBarycenter = [0, 0];
    let totalWeight = 0;
    
    if (weighted) {
        // Use edge-length-weighted barycenter (as in paper)
        const { edgeLengths, edgeMidpoints } = calculateEdgeProperties(projectedVertices, edges);
        
        for (let i = 0; i < edges.length; i++) {
            const weight = edgeLengths[i];
            currentBarycenter[0] += edgeMidpoints[i][0] * weight;
            currentBarycenter[1] += edgeMidpoints[i][1] * weight;
            totalWeight += weight;
        }
    } else {
        // Use simple vertex average
        for (const v of projectedVertices) {
            currentBarycenter[0] += v[0];
            currentBarycenter[1] += v[1];
        }
        totalWeight = projectedVertices.length;
    }
    
    if (totalWeight > 0) {
        currentBarycenter[0] /= totalWeight;
        currentBarycenter[1] /= totalWeight;
    }
    
    // Calculate translation needed
    const dx = targetBarycenter[0] - currentBarycenter[0];
    const dy = targetBarycenter[1] - currentBarycenter[1];
    
    // Apply translation to all vertices
    for (let i = 0; i < projectedVertices.length; i++) {
        projectedVertices[i][0] += dx;
        projectedVertices[i][1] += dy;
    }
    
    return projectedVertices;
}

/**
 * Applies total length constraint by uniformly scaling the curve
 * @param {Array} vertices - Vertex positions 
 * @param {Array} edges - Edge connections
 * @param {Number} targetLength - Target total curve length
 * @returns {Array} - Vertices with scaled total length
 */
export function enforceTotalLength(vertices, edges, targetLength) {
    const projectedVertices = vertices.map(v => [...v]);
    
    // Calculate curve centroid to scale from
    const centroid = [0, 0];
    for (const v of projectedVertices) {
        centroid[0] += v[0];
        centroid[1] += v[1];
    }
    centroid[0] /= projectedVertices.length;
    centroid[1] /= projectedVertices.length;
    
    // Calculate current total length
    let currentLength = 0;
    for (const edge of edges) {
        const v1 = projectedVertices[edge[0]];
        const v2 = projectedVertices[edge[1]];
        const dx = v2[0] - v1[0];
        const dy = v2[1] - v1[1];
        currentLength += Math.sqrt(dx*dx + dy*dy);
    }
    
    // Calculate scale factor
    if (currentLength < 1e-8) return projectedVertices; // Avoid division by zero
    const scale = targetLength / currentLength;
    
    // Apply uniform scaling from centroid
    for (let i = 0; i < projectedVertices.length; i++) {
        projectedVertices[i][0] = centroid[0] + (projectedVertices[i][0] - centroid[0]) * scale;
        projectedVertices[i][1] = centroid[1] + (projectedVertices[i][1] - centroid[1]) * scale;
    }
    
    return projectedVertices;
}

/**
 * Applies multiple constraints in sequence
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {Object} constraints - Constraint configuration
 * @param {Object} additionalData - Additional data for constraints
 * @returns {Array} - Vertices after applying all constraints
 */
export function applyConstraints(vertices, edges, constraints, additionalData = {}) {
    const { initialEdgeLengths = null } = additionalData;
    let result = vertices.map(v => [...v]);
    
    // Apply constraints in a specific order for stability
    
    // 1. First, maintain edge lengths (if enabled)
    if (constraints.maintainEdgeLengths && initialEdgeLengths) {
        console.log("Applying edge length preservation");
        result = projectConstraints(result, edges, initialEdgeLengths);
    }
    
    // 2. Then enforce total length (if enabled)
    if (constraints.totalLength) {
        const targetLength = constraints.targetTotalLength || 
            (initialEdgeLengths ? initialEdgeLengths.reduce((a, b) => a + b, 0) : null);
        
        if (targetLength) {
            console.log(`Applying total length constraint: ${targetLength}`);
            result = enforceTotalLength(result, edges, targetLength);
        }
    }
    
    // 3. Finally, enforce barycenter (if enabled)
    if (constraints.barycenter) {
        console.log(`Applying barycenter constraint: ${constraints.barycenterTarget}`);
        result = enforceBarycenter(result, edges, {
            targetBarycenter: constraints.barycenterTarget || [0, 0],
            weighted: true
        });
    }
    
    return result;
}

/**
 * Projects gradient onto the tangent space of active constraints
 * @param {Array} gradient - Gradient to project
 * @param {Array} vertices - Current vertex positions
 * @param {Array} edges - Edge connections
 * @param {Object} constraints - Active constraints configuration
 * @returns {Array} - Projected gradient
 */
export function projectGradientOntoConstraints(gradient, vertices, edges, constraints) {
    // For simplicity, we'll just zero out components
    // that would violate important constraints
    
    let projectedGradient = gradient.map(g => [...g]);
    
    // If barycenter constraint is active, remove translation component
    if (constraints.barycenter) {
        // Calculate average gradient (translation component)
        const avgGradient = [0, 0];
        for (const g of projectedGradient) {
            avgGradient[0] += g[0] / projectedGradient.length;
            avgGradient[1] += g[1] / projectedGradient.length;
        }
        
        // Subtract translation component from all gradients
        for (let i = 0; i < projectedGradient.length; i++) {
            projectedGradient[i][0] -= avgGradient[0];
            projectedGradient[i][1] -= avgGradient[1];
        }
    }
    
    return projectedGradient;
}