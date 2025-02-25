// src/lib/optimization.js
import {
    calculateDifferential,
    calculateEdgeProperties,
    calculateDiscreteEnergy
} from '$lib/energyCalculations';
import { computePreconditionedGradient } from '$lib/innerProduct';
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from '$lib/stores';

// Configuration toggles
const usePreconditioned = false; // Use preconditioned gradient by default
const applyProjectConstraints = false; // Enable length constraints
const applyBarycenter = true; // Enable barycenter constraint

function l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths) {
    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    const gradient = differential.map(([dx, dy]) => [dx, dy]);
    
    // Normalize the gradient to avoid too large steps
    // const gradNorm = Math.sqrt(gradient.flat().reduce((sum, val) => sum + val * val, 0)) || 1;
    const stepSize = get(config).l2StepSize; // Get step size from config
    
    // Create new vertices by taking a step in the negative gradient direction
    const newVertices = vertices.map((vertex, i) => [
        vertex[0] - stepSize * gradient[i][0] ,
        vertex[1] - stepSize * gradient[i][1]
    ]);

    // Apply constraints if enabled
    let constrainedVertices = [...newVertices];
    if (applyProjectConstraints) {
        constrainedVertices = projectConstraints(constrainedVertices, edges, initialEdgeLengths);
    }
    
    
    if (applyBarycenter) {
        enforceBarycenter(constrainedVertices);
    }

    return constrainedVertices;
}

function preconditionedGradientDescentStep(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    initialEdgeLengths
) {
    const { edgeTangents } = calculateEdgeProperties(vertices, edges);
    
    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    
    // Compute the preconditioned gradient
    let gradient;
    try {
        gradient = computePreconditionedGradient(
            vertices,
            edges,
            edgeTangents,
            alpha,
            beta,
            differential
        );
    } catch (e) {
        console.warn('Preconditioned gradient failed, falling back to L2 gradient:', e);
        return l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths);
    }




    const d = gradient.map(([gx, gy]) => [-gx, -gy]);
    const d_norm = Math.sqrt(d.flat().reduce((sum, val) => sum + val * val, 0)) || 1;
    const d_normalized = d.map(([dx, dy]) => [dx / d_norm, dy / d_norm]);

    const differentialFlat = differential.flat();
    const dFlat = d_normalized.flat();
    const slope = math.dot(differentialFlat, dFlat);
    
    const E_old = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    
    // Line search parameters
    const a_const = get(config).aConst;
    const b_const = get(config).bConst;
    const max_line_search = get(config).maxLineSearch;
    let t = get(config).tauInitial;
    const stepSize = get(config).precondStepSize; // Get step size from config

    // Perform line search to find a good step size
    let vertices_new = [...vertices];
    for (let i = 0; i < max_line_search; i++) {
        vertices_new = vertices.map((vertex, idx) => [
            vertex[0] + t * d_normalized[idx][0] * stepSize,
            vertex[1] + t * d_normalized[idx][1] * stepSize
        ]);
        
        // Apply constraints
        if (applyProjectConstraints) {
            vertices_new = projectConstraints(vertices_new, edges, initialEdgeLengths);
        }
        
        if (applyBarycenter) {
            enforceBarycenter(vertices_new);
        }
        
        const E_new = calculateDiscreteEnergy(vertices_new, edges, alpha, beta, disjointPairs);
        console.log(`Line search iteration ${i}: t=${t}, E_new=${E_new}, E_old=${E_old}, condition=${E_old + a_const * t * slope}`);
        if (E_new <= E_old + a_const * t * slope) {
            console.log('Line search converged at t=', t, 'with energy reduction:', E_old - E_new);
            return vertices_new;
        }
        
        // Backtracking: reduce step size
        t *= b_const;
    }

    console.warn('Line search did not converge, using smallest step size');
    return vertices_new;
}

// Improved constraint projection to maintain edge lengths
function projectConstraints(
    vertices,
    edges,
    initialEdgeLengths,
    maxIterations = get(config).maxConstraintIterations,
    tolerance = get(config).constraintTolerance
) {
    // Create a deep copy to avoid modifying the input
    const projectedVertices = vertices.map(v => [...v]);
    
    for (let iter = 0; iter < maxIterations; iter++) {
        let maxError = 0;
        
        for (let i = 0; i < edges.length; i++) {
            const [v1Idx, v2Idx] = edges[i];
            const v1 = projectedVertices[v1Idx];
            const v2 = projectedVertices[v2Idx];
            
            // Calculate current edge vector and length
            const dx = v2[0] - v1[0];
            const dy = v2[1] - v1[1];
            const currentLength = Math.sqrt(dx * dx + dy * dy) + get(config).epsilonStability;

            // Calculate error relative to target length
            const targetLength = initialEdgeLengths[i];
            const error = (currentLength - targetLength) / targetLength;
            
            // Update the maximum error
            maxError = Math.max(maxError, Math.abs(error));
            
            // If error is significant, adjust vertex positions
            if (Math.abs(error) > tolerance) {
                const correction = error / 2; // Split the correction between vertices
                const scaleFactor = 1 - correction;
                
                // Calculate midpoint
                const midX = (v1[0] + v2[0]) / 2;
                const midY = (v1[1] + v2[1]) / 2;
                
                // Adjust vertex positions to preserve midpoint
                const halfDx = dx / 2;
                const halfDy = dy / 2;
                
                projectedVertices[v1Idx][0] = midX - halfDx * scaleFactor;
                projectedVertices[v1Idx][1] = midY - halfDy * scaleFactor;
                projectedVertices[v2Idx][0] = midX + halfDx * scaleFactor;
                projectedVertices[v2Idx][1] = midY + halfDy * scaleFactor;
            }
        }
        
        // If all constraints are satisfied to tolerance, exit early
        if (maxError < tolerance) {
            console.log(`Constraint projection converged after ${iter + 1} iterations`);
            break;
        }
    }
    
    return projectedVertices;
}

// Improved barycenter constraint
// Improved barycenter constraint with options
function enforceBarycenter(vertices, options = {}) {
    // Default options
    const {
        targetBarycenter = null,  // If null, maintain original position
        centerAtOrigin = false    // If true, center at (0,0) (original behavior)
    } = options;
    
    // Calculate the current center of mass
    const currentBarycenter = [0, 0];
    const n = vertices.length;
    
    for (const vertex of vertices) {
        currentBarycenter[0] += vertex[0] / n;
        currentBarycenter[1] += vertex[1] / n;
    }
    
    // Determine target position based on options
    let target;
    if (centerAtOrigin) {
        target = [0, 0]; // Original behavior (centers at origin)
    } else if (targetBarycenter) {
        target = targetBarycenter; // Use specified target
    } else {
        target = currentBarycenter; // Maintain current position (no movement)
    }
    
    // Calculate the translation needed
    const dx = target[0] - currentBarycenter[0];
    const dy = target[1] - currentBarycenter[1];
    
    // Apply the translation to all vertices
    for (const vertex of vertices) {
        vertex[0] += dx;
        vertex[1] += dy;
    }
    
    return vertices;
}
export function gradientDescentStep(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    initialEdgeLengths
) {
    return usePreconditioned
        ? preconditionedGradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths)
        : l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths);
}

export function createOptimizer(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    maxIterations,
    onUpdate,
    initialEdgeLengths
) {
    if (typeof onUpdate !== 'function') throw new Error('onUpdate must be a function');

    let currentIteration = 0;
    let intervalId = null;
    let lastEnergy = null;
    let stuckCounter = 0;
    
    // Only initialize energy tracking if perturbation is enabled
    const applyPerturbation = get(config).applyPerturbation;
    if (applyPerturbation) {
        lastEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    }

    const optimizer = {
        step: () => {
            if (currentIteration < maxIterations) {
                const newVertices = gradientDescentStep(
                    vertices,
                    edges,
                    alpha,
                    beta,
                    disjointPairs,
                    initialEdgeLengths
                );
                vertices.forEach((v, i) => {
                    v[0] = newVertices[i][0];
                    v[1] = newVertices[i][1];
                });
                
                // Only perform energy calculations and stuck detection if perturbation is enabled
                if (applyPerturbation) {
                    // Calculate new energy
                    const newEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
                    const energyChange = newEnergy - lastEnergy;
                    
                    // Check if we're stuck (minimal energy change)
                    if (Math.abs(energyChange) < get(config).minEnergyChange) {
                        stuckCounter++;
                        
                        // Apply perturbation if stuck for too many iterations
                        if (stuckCounter > get(config).maxStuckIterations) {
                            console.log('Optimizer stuck, applying random perturbation');
                            applyRandomPerturbation(vertices, get(config).perturbationScale);
                            stuckCounter = 0;
                        }
                    } else {
                        stuckCounter = 0;
                    }
                    
                    lastEnergy = newEnergy;
                }
                
                currentIteration++;
                onUpdate();
            } else {
                optimizer.stop();
            }
        },
        start: () => {
            currentIteration = 0;
            stuckCounter = 0;
            
            // Initialize energy tracking on start if perturbation is enabled
            if (applyPerturbation) {
                lastEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
            }
            
            if (!intervalId) intervalId = setInterval(optimizer.step, 20);
        },
        stop: () => {
            if (intervalId) {
                clearInterval(intervalId);
                intervalId = null;
            }
        }
    };

    return optimizer;
}

// Helper function to apply a small random perturbation when optimization gets stuck
function applyRandomPerturbation(vertices, scale) {
    for (const vertex of vertices) {
        vertex[0] += (Math.random() - 0.5) * scale;
        vertex[1] += (Math.random() - 0.5) * scale;
    }
}