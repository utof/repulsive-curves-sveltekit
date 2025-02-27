// src/lib/optimization.js
import {
    calculateDifferential,
    calculateEdgeProperties,
    calculateDiscreteEnergy,
    calculateDiscreteEnergyWithSubvertices
} from '$lib/energyCalculations';
import { computePreconditionedGradient } from '$lib/innerProduct';
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from '$lib/stores';
import { updateKernelState } from '$lib/graphState'; // Import to update subvertices

// Import only the constraints we need
import { projectConstraints, enforceBarycenter } from '$lib/constraints';

// Configuration toggles
const usePreconditioned = true;
const useEdgeLengthConstraint = false;  // Only apply edge length constraint
const useBarycenterConstraint = false; // Don't use barycenter constraint

function l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths) {
    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    const gradient = differential.map(([dx, dy]) => [dx, dy]);
    
    const stepSize = get(config).l2StepSize;
    
    const newVertices = vertices.map((vertex, i) => [
        vertex[0] - stepSize * gradient[i][0],
        vertex[1] - stepSize * gradient[i][1]
    ]);

    let constrainedVertices = [...newVertices];
    if (useEdgeLengthConstraint) {
        constrainedVertices = projectConstraints(constrainedVertices, edges, initialEdgeLengths);
    }
    
    if (useBarycenterConstraint) {
        constrainedVertices = enforceBarycenter(constrainedVertices);
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

    const stepSize = get(config).precondStepSize;
    const useLineSearch = get(config).useLineSearch;

    // If line search is disabled, just take a fixed step
    if (!useLineSearch) {
        const newVertices = vertices.map((vertex, idx) => [
            vertex[0] + d[idx][0] * stepSize,
            vertex[1] + d[idx][1] * stepSize
        ]);
        
        let constrainedVertices = [...newVertices];
        if (useEdgeLengthConstraint) {
            constrainedVertices = projectConstraints(constrainedVertices, edges, initialEdgeLengths);
        }
        
        if (useBarycenterConstraint) {
            constrainedVertices = enforceBarycenter(constrainedVertices);
        }
        
        return constrainedVertices;
    }

    // Line search is enabled, so we perform backtracking line search
    const differentialFlat = differential.flat();
    const dFlat = d.flat();
    const slope = math.dot(differentialFlat, dFlat);
    
    const E_old = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    
    const a_const = get(config).aConst;
    const b_const = get(config).bConst;
    const max_line_search = get(config).maxLineSearch;
    let t = get(config).tauInitial;

    let vertices_new = [...vertices];
    for (let i = 0; i < max_line_search; i++) {
        vertices_new = vertices.map((vertex, idx) => [
            vertex[0] + t * d[idx][0] * stepSize,
            vertex[1] + t * d[idx][1] * stepSize
        ]);
        
        if (useEdgeLengthConstraint) {
            vertices_new = projectConstraints(vertices_new, edges, initialEdgeLengths);
        }
        
        if (useBarycenterConstraint) {
            vertices_new = enforceBarycenter(vertices_new);
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