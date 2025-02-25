// src/lib/energyCalculations.js
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from '$lib/stores';

let logging = false;

export function calculateEdgeProperties(vertices, edges) {
    const dim = get(config).dim;
    const edgeLengths = [];
    const edgeTangents = [];
    const edgeMidpoints = [];

    for (const edge of edges) {
        const v1 = vertices[edge[0]];
        const v2 = vertices[edge[1]];

        const diff = v2.map((coord, i) => coord - v1[i]);
        const length = Math.sqrt(diff.reduce((sum, d) => sum + d * d, 0));
        edgeLengths.push(length);

        const unitTangent = length > 0 ? diff.map(d => d / length) : new Array(dim).fill(0);
        edgeTangents.push(unitTangent);

        const midpoint = v1.map((coord, i) => (coord + v2[i]) / 2);
        edgeMidpoints.push(midpoint);

        if (logging) {
            console.log(`Edge [${edge[0]}, ${edge[1]}]: length = ${length}, tangent = ${unitTangent}, midpoint = ${midpoint}`);
        }
    }

    return { edgeLengths, edgeTangents, edgeMidpoints };
}

export function tangentPointKernel(p, q, T, alpha, beta) {
    const dim = get(config).dim;
    const epsilon = get(config).epsilonKernel;

    const diff = p.map((coord, i) => coord - q[i]);
    const diffNorm = Math.sqrt(diff.reduce((sum, d) => sum + d * d, 0)) + epsilon;

    let cross;
    if (dim === 2) {
        cross = Math.abs(T[0] * diff[1] - T[1] * diff[0]); // 2D determinant
    } else if (dim === 3) {
        const crossVec = [
            T[1] * diff[2] - T[2] * diff[1],
            T[2] * diff[0] - T[0] * diff[2],
            T[0] * diff[1] - T[1] * diff[0]
        ];
        cross = Math.sqrt(crossVec.reduce((sum, c) => sum + c * c, 0)); // 3D cross product norm
    } else {
        throw new Error('Unsupported dimension');
    }

    const numerator = Math.pow(cross, alpha);
    const denominator = Math.pow(diffNorm, beta);
    const result = numerator / denominator;

    if (!isFinite(result)) {
        console.warn('Invalid kernel result:', result, 'from inputs:', p, q, T, 'with cross:', cross, 'diffNorm:', diffNorm);
        return 0;
    }

    return result;
}

export function calculateDisjointEdgePairs(edges) {
    const numEdges = edges.length;
    const disjointPairs = [];

    for (let i = 0; i < numEdges; i++) {
        disjointPairs[i] = [];
        for (let j = 0; j < numEdges; j++) {
            if (i === j) continue;

            const edge1 = edges[i];
            const edge2 = edges[j];

            if (
                edge1[0] !== edge2[0] &&
                edge1[0] !== edge2[1] &&
                edge1[1] !== edge2[0] &&
                edge1[1] !== edge2[1]
            ) {
                disjointPairs[i].push(j);
            }
        }
    }
    return disjointPairs;
}

export function calculateDiscreteKernel(vertices, edges, edgeTangents, alpha, beta, disjointPairs) {
    const numEdges = edges.length;
    const kernelMatrix = math.zeros(numEdges, numEdges);

    if (!disjointPairs || !Array.isArray(disjointPairs) || disjointPairs.length === 0) {
        console.warn('No disjoint pairs found, returning zero kernel matrix');
        return kernelMatrix;
    }

    for (let i = 0; i < numEdges; i++) {
        if (!disjointPairs[i]) continue;

        for (const j of disjointPairs[i]) {
            if (i < edges.length && j < edges.length) {
                let sum = 0;
                const combinations = [
                    [vertices[edges[i][0]], vertices[edges[j][0]]],
                    [vertices[edges[i][0]], vertices[edges[j][1]]],
                    [vertices[edges[i][1]], vertices[edges[j][0]]],
                    [vertices[edges[i][1]], vertices[edges[j][1]]]
                ];

                for (const [p, q] of combinations) {
                    sum += tangentPointKernel(p, q, edgeTangents[i], alpha, beta);
                }
                kernelMatrix.set([i, j], sum / 4);
                kernelMatrix.set([j, i], sum / 4);
            }
        }
    }
    return kernelMatrix;
}

export function calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs) {
    const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
    const kernelMatrix = calculateDiscreteKernel(vertices, edges, edgeTangents, alpha, beta, disjointPairs);

    let totalEnergy = 0;
    const numEdges = edges.length;

    for (let i = 0; i < numEdges; i++) {
        for (const j of disjointPairs[i]) {
            if (i < edges.length && j < edges.length) {
                const kernelValue = kernelMatrix.get([i, j]);
                totalEnergy += kernelValue * edgeLengths[i] * edgeLengths[j];
            }
        }
    }
    return totalEnergy / 2;
}

export function calculateDifferential(vertices, edges, alpha, beta, disjointPairs) {
    const dim = get(config).dim;
    const method = get(config).differentialMethod;

    if (method === 'finiteDifference') {
        const h = get(config).finiteDiffH;
        const numVertices = vertices.length;
        const differential = [];

        const E_original = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

        for (let i = 0; i < numVertices; i++) {
            differential[i] = new Array(dim).fill(0);
            for (let d = 0; d < dim; d++) {
                const vertices_perturbed = vertices.map(v => [...v]);
                vertices_perturbed[i][d] += h;
                const E_perturbed = calculateDiscreteEnergy(vertices_perturbed, edges, alpha, beta, disjointPairs);
                differential[i][d] = (E_perturbed - E_original) / h;
            }
        }
        return differential;
    } else {
        throw new Error('Analytical differential not implemented for general dimensions');
    }
}