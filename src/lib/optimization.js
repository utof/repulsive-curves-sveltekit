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
const usePreconditioned = 1; // 0 for L2, 1 for Preconditioned
const applyProjectConstraints = false; // Comment out to disable projectConstraints
const applyBarycenter = true;         // Comment out to disable enforceBarycenter

function l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths) {
    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    const gradient = differential.map(([dx, dy]) => [dx, dy]);
    const stepSize = 1000;
    const newVertices = vertices.map((vertex, i) => [
        vertex[0] - stepSize * gradient[i][0],
        vertex[1] - stepSize * gradient[i][1]
    ]);

    // Apply constraints based on toggles
    if (applyProjectConstraints) projectConstraints(newVertices, edges, initialEdgeLengths);
    if (applyBarycenter) enforceBarycenter(newVertices, edges, initialEdgeLengths);

    return newVertices;
}

function preconditionedGradientDescentStep(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    initialEdgeLengths,
    a_const = get(config).aConst,
    b_const = get(config).bConst,
    max_line_search = get(config).maxLineSearch
) {
    const stepsize = 10;
    const { edgeTangents, edgeLengths } = calculateEdgeProperties(vertices, edges);
    console.log('Gradient Descent Step - Initial vertices:', JSON.stringify(vertices));

    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    const gradient = computePreconditionedGradient(
        vertices,
        edges,
        edgeTangents,
        alpha,
        beta,
        differential
    );

    const d = gradient.map(([gx, gy]) => [-gx, -gy]);
    const d_norm = Math.sqrt(d.flat().reduce((sum, val) => sum + val * val, 0)) || 1;
    const d_normalized = d.map(([dx, dy]) => [dx / d_norm, dy / d_norm]);

    const differentialFlat = differential.flat();
    const dFlat = d_normalized.flat();
    const slope = math.dot(differentialFlat, dFlat);

    const E_old = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

    let t = get(config).tauInitial;
    for (let i = 0; i < max_line_search; i++) {
        const vertices_new = vertices.map((vertex, idx) => [
            vertex[0] + t * d_normalized[idx][0] * stepsize,
            vertex[1] + t * d_normalized[idx][1] * stepsize
        ]);

        // Apply constraints based on toggles
        if (applyProjectConstraints) projectConstraints(vertices_new, edges, initialEdgeLengths);
        if (applyBarycenter) enforceBarycenter(vertices_new, edges, edgeLengths);

        const E_new = calculateDiscreteEnergy(vertices_new, edges, alpha, beta, disjointPairs);
        console.log(`Line search iteration ${i}: t=${t}, E_new=${E_new}, E_old=${E_old}, condition=${E_old + a_const * t * slope}`);
        if (E_new <= E_old + a_const * t * slope) {
            console.log('Line search converged at t=', t);
            return vertices_new;
        }
        // t *= b_const;
    }

    console.warn('Line search did not converge after max iterations');
    return vertices;
}

function projectConstraints(
    vertices,
    edges,
    initialEdgeLengths,
    maxIterations = get(config).maxLineSearch,
    tolerance = get(config).constraintTolerance
) {
    for (let iter = 0; iter < maxIterations; iter++) {
        let maxError = 0;
        for (let I = 0; I < edges.length; I++) {
            const [i, j] = edges[I];
            const vi = vertices[i];
            const vj = vertices[j];
            const d = [vj[0] - vi[0], vj[1] - vi[1]];
            const currentLength = Math.sqrt(d[0] * d[0] + d[1] * d[1]);
            const targetLength = initialEdgeLengths[I];
            const error = currentLength - targetLength;
            if (Math.abs(error) > tolerance) {
                const correction = (error / currentLength) * 0.5;
                const delta = [d[0] * correction, d[1] * correction];
                vertices[i][0] += delta[0];
                vertices[i][1] += delta[1];
                vertices[j][0] -= delta[0];
                vertices[j][1] -= delta[1];
                maxError = Math.max(maxError, Math.abs(error));
            }
        }
        if (maxError < tolerance) break;
    }
}

function enforceBarycenter(vertices, edges, edgeLengths) {
    const totalLength = edgeLengths.reduce((sum, l) => sum + l, 0);
    let barycenter = [0, 0];
    for (let i = 0; i < edges.length; i++) {
        const [v0, v1] = edges[i];
        const midpoint = [
            (vertices[v0][0] + vertices[v1][0]) / 2,
            (vertices[v0][1] + vertices[v1][1]) / 2
        ];
        barycenter[0] += (edgeLengths[i] / totalLength) * midpoint[0];
        barycenter[1] += (edgeLengths[i] / totalLength) * midpoint[1];
    }
    vertices.forEach((v) => {
        v[0] -= barycenter[0];
        v[1] -= barycenter[1];
    });
}

export function gradientDescentStep(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    initialEdgeLengths,
    a_const = get(config).aConst,
    b_const = get(config).bConst,
    max_line_search = get(config).maxLineSearch
) {
    return usePreconditioned
        ? preconditionedGradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths, a_const, b_const, max_line_search)
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
                currentIteration++;
                onUpdate();
            } else {
                optimizer.stop();
            }
        },
        start: () => {
            currentIteration = 0;
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