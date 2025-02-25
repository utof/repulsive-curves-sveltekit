// src/lib/optimization.js
import { calculateDifferential, calculateEdgeProperties, calculateDiscreteEnergy } from '$lib/energyCalculations';
import { computePreconditionedGradient } from '$lib/innerProduct';
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from '$lib/stores';
import { updateKernelState } from '$lib/graphState';

const usePreconditioned = true;
const applyProjectConstraints = false;
const applyBarycenter = false;

function l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths) {
    const dim = get(config).dim;
    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    const gradient = differential.map(coord => [...coord]);

    const stepSize = get(config).l2StepSize;

    const newVertices = vertices.map((vertex, i) =>
        vertex.map((coord, d) => coord - stepSize * gradient[i][d])
    );

    let constrainedVertices = [...newVertices];
    if (applyProjectConstraints) {
        constrainedVertices = projectConstraints(constrainedVertices, edges, initialEdgeLengths);
    }

    if (applyBarycenter) {
        enforceBarycenter(constrainedVertices);
    }

    return constrainedVertices;
}

function preconditionedGradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths) {
    const dim = get(config).dim;
    const { edgeTangents } = calculateEdgeProperties(vertices, edges);

    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);

    let gradient;
    try {
        gradient = computePreconditionedGradient(vertices, edges, edgeTangents, alpha, beta, differential);
    } catch (e) {
        console.warn('Preconditioned gradient failed, falling back to L2 gradient:', e);
        return l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths);
    }

    const d = gradient.map(coord => coord.map(val => -val));
    const d_norm = Math.sqrt(d.flat().reduce((sum, val) => sum + val * val, 0)) || 1;
    const d_normalized = d.map(coord => coord.map(val => val / d_norm));

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
        vertices_new = vertices.map((vertex, idx) =>
            vertex.map((coord, d) => coord + t * d_normalized[idx][d] * stepSize)
        );

        if (applyProjectConstraints) {
            vertices_new = projectConstraints(vertices_new, edges, initialEdgeLengths);
        }

        if (applyBarycenter) {
            enforceBarycenter(vertices_new);
        }

        const E_new = calculateDiscreteEnergy(vertices_new, edges, alpha, beta, disjointPairs);
        if (E_new <= E_old + a_const * t * slope) {
            return vertices_new;
        }

        t *= b_const;
    }

    console.warn('Line search did not converge, using smallest step size');
    return vertices_new;
}

function projectConstraints(vertices, edges, initialEdgeLengths) {
    const dim = get(config).dim;
    const projectedVertices = vertices.map(v => [...v]);
    const maxIterations = get(config).maxConstraintIterations;
    const tolerance = get(config).constraintTolerance;

    for (let iter = 0; iter < maxIterations; iter++) {
        let maxError = 0;

        for (let i = 0; i < edges.length; i++) {
            const [v1Idx, v2Idx] = edges[i];
            const v1 = projectedVertices[v1Idx];
            const v2 = projectedVertices[v2Idx];

            const diff = v2.map((coord, d) => coord - v1[d]);
            const currentLength = Math.sqrt(diff.reduce((sum, d) => sum + d * d, 0)) + get(config).epsilonStability;

            const targetLength = initialEdgeLengths[i];
            const error = (currentLength - targetLength) / targetLength;

            maxError = Math.max(maxError, Math.abs(error));

            if (Math.abs(error) > tolerance) {
                const correction = error / 2;
                const scaleFactor = 1 - correction;

                const mid = v1.map((coord, d) => (coord + v2[d]) / 2);
                const halfDiff = diff.map(d => d / 2);

                projectedVertices[v1Idx] = mid.map((m, d) => m - halfDiff[d] * scaleFactor);
                projectedVertices[v2Idx] = mid.map((m, d) => m + halfDiff[d] * scaleFactor);
            }
        }

        if (maxError < tolerance) {
            break;
        }
    }

    return projectedVertices;
}

function enforceBarycenter(vertices, options = {}) {
    const dim = get(config).dim;
    const { targetBarycenter = null, centerAtOrigin = false } = options;

    const currentBarycenter = new Array(dim).fill(0);
    const n = vertices.length;

    for (const vertex of vertices) {
        for (let d = 0; d < dim; d++) {
            currentBarycenter[d] += vertex[d] / n;
        }
    }

    let target = targetBarycenter || (centerAtOrigin ? new Array(dim).fill(0) : currentBarycenter);
    const delta = target.map((t, d) => t - currentBarycenter[d]);

    for (const vertex of vertices) {
        for (let d = 0; d < dim; d++) {
            vertex[d] += delta[d];
        }
    }

    return vertices;
}

export function gradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths) {
    return usePreconditioned
        ? preconditionedGradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths)
        : l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths);
}

export function createOptimizer(vertices, edges, alpha, beta, disjointPairs, maxIterations, onUpdate, initialEdgeLengths) {
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
                const newVertices = gradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths);
                vertices.forEach((v, i) => v.splice(0, v.length, ...newVertices[i]));

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
    const dim = get(config).dim;
    for (const vertex of vertices) {
        for (let d = 0; d < dim; d++) {
            vertex[d] += (Math.random() - 0.5) * scale;
        }
    }
}