// src/lib/innerProduct.js
import * as math from 'mathjs';
import { calculateEdgeProperties, calculateDisjointEdgePairs } from '$lib/energyCalculations';
import { get } from 'svelte/store';
import { config } from '$lib/stores';

function build_weights(alpha, beta, edges, disjointPairs, vertices, edgeTangents, edgeLengths) {
    const dim = get(config).dim;
    const edge_num = edges.length;
    const s = (beta - 1) / alpha;
    const sigma = s - 1;

    const W = math.zeros(edge_num, edge_num);
    const W0 = math.zeros(edge_num, edge_num);

    for (let I = 0; I < disjointPairs.length; I++) {
        for (const J of disjointPairs[I]) {
            let elt1 = 0;
            let elt2 = 0;

            for (let a = 0; a < 2; a++) {
                for (let b = 0; b < 2; b++) {
                    const i = edges[I][a];
                    const j = edges[J][b];
                    const p = vertices[i];
                    const q = vertices[j];
                    const diff = p.map((coord, k) => coord - q[k]);
                    const epsilon = get(config).epsilonStability;
                    let diff_norm = Math.sqrt(diff.reduce((sum, d) => sum + d * d, 0)) + epsilon;

                    const term1 = 1 / Math.pow(diff_norm, 2 * sigma + 1);
                    elt1 += term1;

                    const alph = 2;
                    const bet = 4;

                    let cross_norm;
                    if (dim === 2) {
                        const cross = diff[0] * edgeTangents[I][1] - diff[1] * edgeTangents[I][0];
                        cross_norm = Math.abs(cross);
                    } else if (dim === 3) {
                        const crossVec = [
                            diff[1] * edgeTangents[I][2] - diff[2] * edgeTangents[I][1],
                            diff[2] * edgeTangents[I][0] - diff[0] * edgeTangents[I][2],
                            diff[0] * edgeTangents[I][1] - diff[1] * edgeTangents[I][0]
                        ];
                        cross_norm = Math.sqrt(crossVec.reduce((sum, c) => sum + c * c, 0));
                    }

                    const k_numerator = Math.pow(cross_norm, alph);
                    const k_denominator = Math.pow(diff_norm, bet);
                    const k = k_numerator / k_denominator;

                    const term2 = k / Math.pow(diff_norm, 2 * sigma + 1);
                    elt2 += term2;
                }
            }

            const w_ij_factor = 0.25 * edgeLengths[I] * edgeLengths[J];
            W.set([I, J], w_ij_factor * elt1);
            W0.set([I, J], w_ij_factor * elt2);
        }
    }

    return { W, W0 };
}

function calculateLowOrderTerm(vertices, edges, W0) {
    const numVertices = vertices.length;
    const B0 = math.zeros(numVertices, numVertices);
    const disjointEdges = calculateDisjointEdgePairs(edges);

    for (let I = 0; I < edges.length; I++) {
        for (const J of disjointEdges[I]) {
            const w_IJ_0 = W0.get([I, J]);

            for (let a = 0; a < 2; a++) {
                for (let b = 0; b < 2; b++) {
                    const i_a = edges[I][a];
                    const i_b = edges[I][b];
                    const j_a = edges[J][a];
                    const j_b = edges[J][b];

                    B0.set([i_a, i_b], B0.get([i_a, i_b]) + 0.25 * w_IJ_0);
                    B0.set([j_a, j_b], B0.get([j_a, j_b]) + 0.25 * w_IJ_0);
                    B0.set([i_a, j_b], B0.get([i_a, j_b]) - 0.25 * w_IJ_0);
                    B0.set([j_a, i_b], B0.get([j_a, i_b]) - 0.25 * w_IJ_0);
                }
            }
        }
    }
    return B0;
}

function calculateHighOrderTerm(vertices, edges, W, edgeLengths, edgeTangents) {
    const dim = get(config).dim;
    const numVertices = vertices.length;
    const B = math.zeros(numVertices, numVertices);
    const disjointEdges = calculateDisjointEdgePairs(edges);

    for (let I = 0; I < edges.length; I++) {
        for (const J of disjointEdges[I]) {
            const l_I = edgeLengths[I];
            const l_J = edgeLengths[J];
            const T_I = edgeTangents[I];
            const T_J = edgeTangents[J];
            const w_IJ = W.get([I, J]);
            const dot_TI_TJ = T_I.reduce((sum, val, idx) => sum + val * T_J[idx], 0);

            for (let a = 0; a < 2; a++) {
                for (let b = 0; b < 2; b++) {
                    const sign = Math.pow(-1, a + b);
                    const i_a = edges[I][a];
                    const i_b = edges[I][b];
                    const j_a = edges[J][a];
                    const j_b = edges[J][b];

                    const val_1 = (sign * w_IJ) / (l_I * l_I);
                    const val_2 = (sign * w_IJ) / (l_J * l_J);
                    const val_3 = (sign * w_IJ * dot_TI_TJ) / (l_I * l_J);

                    B.set([i_a, i_b], B.get([i_a, i_b]) + val_1);
                    B.set([j_a, j_b], B.get([j_a, j_b]) + val_2);
                    B.set([i_a, j_b], B.get([i_a, j_b]) - val_3);
                    B.set([j_a, i_b], B.get([j_a, i_b]) - val_3);
                }
            }
        }
    }
    return B;
}

export function build_A_bar(alpha, beta, vertices, edges) {
    const dim = get(config).dim;
    const numVertices = vertices.length;
    const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
    const disjointEdges = calculateDisjointEdgePairs(edges);

    const { W, W0 } = build_weights(alpha, beta, edges, disjointEdges, vertices, edgeTangents, edgeLengths);
    const B = calculateHighOrderTerm(vertices, edges, W, edgeLengths, edgeTangents);
    const B0 = calculateLowOrderTerm(vertices, edges, W0);
    let A = math.add(B, B0);

    const epsilon = get(config).epsilonStability;
    const reg = math.multiply(epsilon, math.identity(numVertices));
    A = math.add(A, reg);

    const A_bar = math.zeros(dim * numVertices, dim * numVertices);
    for (let d = 0; d < dim; d++) {
        const start = d * numVertices;
        const end = (d + 1) * numVertices;
        A_bar.subset(math.index(math.range(start, end), math.range(start, end)), A);
    }

    return { A_bar, B, B0 };
}

export function computePreconditionedGradient(vertices, edges, edgeTangents, alpha, beta, differential) {
    const dim = get(config).dim;
    const numVertices = vertices.length;
    const { A_bar } = build_A_bar(alpha, beta, vertices, edges);
    const differentialFlat = differential.flat();

    let gradFlat;
    try {
        gradFlat = math.lusolve(A_bar, differentialFlat);
    } catch (e) {
        console.error('Linear solve failed:', e);
        throw new Error('Failed to compute preconditioned gradient due to singular matrix');
    }

    const gradArray = gradFlat.toArray();
    const grad = [];
    for (let i = 0; i < numVertices; i++) {
        grad[i] = gradArray.slice(i * dim, (i + 1) * dim);
    }
    return grad;
}