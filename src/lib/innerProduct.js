// src/lib/innerProduct.js
import * as math from 'mathjs';
import {
	calculateEdgeProperties,
	calculateDisjointEdgePairs,
	tangentPointKernel
} from '$lib/energyCalculations';

/**
 * Builds the weight matrices W and W0 used in energy calculations.
 * @param {number} alpha - Energy parameter.
 * @param {number} beta - Energy parameter.
 * @param {Array<Array<number>>} edges - Array of edges, where each edge is an array of two vertex indices.
 * @param {Array<Array<number>>} disjointPairs - Array of disjoint edge pairs.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edgeTangents - Array of edge tangent vectors.
 * @param {Array<number>} edgeLengths - Array of edge lengths.
 * @returns {{W: math.Matrix, W0: math.Matrix}} - The weight matrices W and W0.
 */
function build_weights(alpha, beta, edges, disjointPairs, vertices, edgeTangents, edgeLengths) {
	const edge_num = edges.length;
	const sigma = (beta - 1) / alpha - 1;

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
					const diff = math.subtract(p, q);
					let diff_norm = math.norm(diff);

					// Handle near-zero distances
					diff_norm = Math.max(diff_norm, 1e-8);

					const term1 = 1 / Math.pow(diff_norm, 2 * sigma + 1);
					elt1 += term1;

					// Constants for B0 calculation
					const alph = 2;
					const bet = 4;

					const cross = math.det([diff, edgeTangents[I]]);
					const cross_norm = math.norm(cross);

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

/**
 * Calculates the low-order term matrix B0.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edges - Array of edges.
 * @param {math.Matrix} W0 - The weight matrix W0.
 * @returns {math.Matrix} - The low-order term matrix B0.
 */
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
	console.log('Low order term B0:', B0.toArray());
	return B0;
}

/**
 * Calculates the high-order term matrix B.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edges - Array of edges.
 * @param {math.Matrix} W - The weight matrix W.
 * @param {Array<number>} edgeLengths - Array of edge lengths.
 * @param {Array<Array<number>>} edgeTangents - Array of edge tangent vectors.
 * @returns {math.Matrix} - The high-order term matrix B.
 */
function calculateHighOrderTerm(vertices, edges, W, edgeLengths, edgeTangents) {
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
			const dot_TI_TJ = math.dot(T_I, T_J);

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

	console.log('High order term B:', B.toArray());
	return B;
}

/**
 * Calculates the discrete inner product matrix A and its components.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edges - Array of edges.
 * @param {number} alpha - Energy parameter.
 * @param {number} beta - Energy parameter.
 * @returns {{A_reg: math.Matrix, B0: math.Matrix, B: math.Matrix}} - The regularized inner product matrix A_reg and component matrices B0 and B.
 */
export function calculateDiscreteInnerProduct(vertices, edges, alpha, beta) {
	console.log('Calculating discrete inner product with alpha:', alpha, 'beta:', beta);
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);

	const { W, W0 } = build_weights(
		alpha,
		beta,
		edges,
		calculateDisjointEdgePairs(edges),
		vertices,
		edgeTangents,
		edgeLengths
	);
	const B = calculateHighOrderTerm(vertices, edges, W, edgeLengths, edgeTangents);
	const B0 = calculateLowOrderTerm(vertices, edges, W0);

	let A = math.add(B0, B);
	console.log('Initial A Matrix (B0 + B):', A.toArray());

	// Enforce constraint: vertex 0 is fixed
	const numVertices = vertices.length;
	for (let j = 0; j < numVertices; j++) {
		A.set([0, j], j === 0 ? 1 : 0);
		A.set([j, 0], j === 0 ? 1 : 0);
	}
	console.log('A Matrix after fixing vertex 0:', A.toArray());

	// Add regularization
	const reg = math.multiply(1e-6, math.identity(numVertices));
	const A_reg = math.add(A, reg);
	console.log('A Matrix with regularization:', A_reg.toArray());

	try {
		const eigenvalues = math.eigs(A_reg).values.toArray();
		console.log('Eigenvalues of A_reg:', eigenvalues);
		const minEigen = Math.min(...eigenvalues);
		console.log('Minimum eigenvalue:', minEigen);
		if (minEigen <= 0) console.warn('A_reg is not positive definite!');
	} catch (e) {
		console.error('Eigenvalue computation failed:', e);
	}

	return { A_reg, B0, B };
}

/**
 * Builds the 2D block-diagonal matrix A_bar.
 * @param {number} alpha - Energy parameter.
 * @param {number} beta - Energy parameter.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edges - Array of edges.
 * @returns {{A_bar: math.Matrix, B: math.Matrix, B0: math.Matrix}} - The block-diagonal matrix A_bar and component matrices B and B0.
 */
export function build_A_bar_2D(alpha, beta, vertices, edges) {
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
	const disjointEdges = calculateDisjointEdgePairs(edges);

	const { W, W0 } = build_weights(
		alpha,
		beta,
		edges,
		disjointEdges,
		vertices,
		edgeTangents,
		edgeLengths
	);
	const B = calculateHighOrderTerm(vertices, edges, W, edgeLengths, edgeTangents);
	const B0 = calculateLowOrderTerm(vertices, edges, W0);
	const A = math.add(B, B0);

	const numVertices = vertices.length;
	const A_bar = math.zeros(2 * numVertices, 2 * numVertices); // 2D, so 2 * numVertices

	// Block-diagonal construction (only two blocks for 2D)
	A_bar.subset(math.index(math.range(0, numVertices), math.range(0, numVertices)), A);
	A_bar.subset(
		math.index(math.range(numVertices, 2 * numVertices), math.range(numVertices, 2 * numVertices)),
		A
	);

	return { A_bar, B, B0 };
}

export function computePreconditionedGradient(
	vertices,
	edges,
	edgeTangents,
	alpha,
	beta,
	differential
) {
	console.log('Computing preconditioned gradient with differential:', differential);
	const numVertices = vertices.length;
	const AFull = math.zeros(numVertices * 2, numVertices * 2);
	const sigma = (beta - 1) / alpha - 1;
	const { A_bar, B, B0 } = build_A_bar_2D(alpha, beta, vertices, edges);
	const A = A_bar;

	console.log('Constructing AFull...');
	for (let i = 0; i < numVertices; i++) {
		for (let j = 0; j < numVertices; j++) {
			const val = A.get([i, j]);
			AFull.set([i * 2, j * 2], val);
			AFull.set([i * 2 + 1, j * 2 + 1], val);
		}
	}
	console.log('AFull Matrix:', AFull.toArray());

	// Flatten the differential
	const differentialFlat = differential.flat();

	let gradFlatFull;
	try {
		// Solve the system A * g = dE  (Equation 20)
		gradFlatFull = math.multiply(math.inv(AFull), differentialFlat);
		console.log('Preconditioned gradient (flat) as matrix:', gradFlatFull.toArray());
	} catch (e) {
		console.error('Inversion failed:', e);
		// Fallback (though ideally we should never reach here)
		gradFlatFull = math.matrix(differentialFlat);
	}

	const gradArray = gradFlatFull.toArray();
	const grad = [];
	for (let i = 0; i < numVertices; i++) {
		grad[i] = [gradArray[i * 2], gradArray[i * 2 + 1]];
		console.log(`Gradient for vertex ${i}: [${grad[i][0]}, ${grad[i][1]}]`);
	}
	return grad;
}
