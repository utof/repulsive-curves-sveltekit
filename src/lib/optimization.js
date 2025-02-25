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
import { updateKernelState } from '$lib/graphState'; // Import to update subvertices

// Configuration toggles
const usePreconditioned = false;
const applyProjectConstraints = false;
const applyBarycenter = true;

function l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths) {
    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    const gradient = differential.map(([dx, dy]) => [dx, dy]);
    
    const stepSize = get(config).l2StepSize;
    
    const newVertices = vertices.map((vertex, i) => [
        vertex[0] - stepSize * gradient[i][0],
        vertex[1] - stepSize * gradient[i][1]
    ]);

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
    
    const a_const = get(config).aConst;
    const b_const = get(config).bConst;
    const max_line_search = get(config).maxLineSearch;
    let t = get(config).tauInitial;
    const stepSize = get(config).precondStepSize;

    let vertices_new = [...vertices];
    for (let i = 0; i < max_line_search; i++) {
        vertices_new = vertices.map((vertex, idx) => [
            vertex[0] + t * d_normalized[idx][0] * stepSize,
            vertex[1] + t * d_normalized[idx][1] * stepSize
        ]);
        
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
        
        t *= b_const;
    }

    console.warn('Line search did not converge, using smallest step size');
    return vertices_new;
}

function projectConstraints(
    vertices,
    edges,
    initialEdgeLengths,
    maxIterations = get(config).maxConstraintIterations,
    tolerance = get(config).constraintTolerance
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
            const currentLength = Math.sqrt(dx * dx + dy * dy) + get(config).epsilonStability;

            const targetLength = initialEdgeLengths[i];
            const error = (currentLength - targetLength) / targetLength;
            
            maxError = Math.max(maxError, Math.abs(error));
            
            if (Math.abs(error) > tolerance) {
                const correction = error / 2;
                const scaleFactor = 1 - correction;
                
                const midX = (v1[0] + v2[0]) / 2;
                const midY = (v1[1] + v2[1]) / 2;
                
                const halfDx = dx / 2;
                const halfDy = dy / 2;
                
                projectedVertices[v1Idx][0] = midX - halfDx * scaleFactor;
                projectedVertices[v1Idx][1] = midY - halfDy * scaleFactor;
                projectedVertices[v2Idx][0] = midX + halfDx * scaleFactor;
                projectedVertices[v2Idx][1] = midY + halfDy * scaleFactor;
            }
        }
        
        if (maxError < tolerance) {
            console.log(`Constraint projection converged after ${iter + 1} iterations`);
            break;
        }
    }
    
    return projectedVertices;
}

function enforceBarycenter(vertices, options = {}) {
    const {
        targetBarycenter = null,
        centerAtOrigin = false
    } = options;
    
    const currentBarycenter = [0, 0];
    const n = vertices.length;
    
    for (const vertex of vertices) {
        currentBarycenter[0] += vertex[0] / n;
        currentBarycenter[1] += vertex[1] / n;
    }
    
    let target;
    if (centerAtOrigin) {
        target = [0, 0];
    } else if (targetBarycenter) {
        target = targetBarycenter;
    } else {
        target = currentBarycenter;
    }
    
    const dx = target[0] - currentBarycenter[0];
    const dy = target[1] - currentBarycenter[1];
    
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
                
                // Update kernel state including subvertices
                updateKernelState(vertices, edges, alpha, beta, disjointPairs);
                
                if (applyPerturbation) {
                    const newEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
                    const energyChange = newEnergy - lastEnergy;
                    
                    if (Math.abs(energyChange) < get(config).minEnergyChange) {
                        stuckCounter++;
                        
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

function applyRandomPerturbation(vertices, scale) {
    for (const vertex of vertices) {
        vertex[0] += (Math.random() - 0.5) * scale;
        vertex[1] += (Math.random() - 0.5) * scale;
    }
}