// src/lib/innerProduct.js
import * as math from 'mathjs';
import {
	calculateEdgeProperties,
	calculateDisjointEdgePairs,
	tangentPointKernel
} from '$lib/energyCalculations';

function build_weights(alpha, beta, edges, disjointPairs, vertices, edgeTangents, edgeLengths) {
	console.log('Entering build_weights function');
	console.log('alpha:', alpha);
	console.log('beta:', beta);
	console.log('edges:', edges);
	console.log('disjointPairs:', disjointPairs);
	console.log('vertices:', vertices);
	console.log('edgeTangents:', edgeTangents);
	console.log('edgeLengths:', edgeLengths);

	const edge_num = edges.length;
	console.log('edge_num:', edge_num);

	const sigma = (beta - 1) / alpha;
	console.log('sigma:', sigma);

	const W = math.zeros(edge_num, edge_num);
	const W0 = math.zeros(edge_num, edge_num);
	console.log('Initialized W:', W);
	console.log('Initialized W0:', W0);

	for (let I = 0; I < disjointPairs.length; I++) {
		console.log('Outer loop: I =', I);
		console.log('disjointPairs[I]:', disjointPairs[I]);

		for (const J_ind of disjointPairs[I]) {
			console.log('  Inner loop: J_ind =', J_ind);
			const J = J_ind;
			console.log('  J =', J);

			let elt1 = 0;
			let elt2 = 0;
			console.log('  Initialized elt1:', elt1);
			console.log('  Initialized elt2:', elt2);

			for (let a = 0; a < 2; a++) {
				console.log('    a-loop: a =', a);
				for (let b = 0; b < 2; b++) {
					console.log('      b-loop: b =', b);

					const i = edges[I][a];
					const j = edges[J][b];
					console.log('      i:', i);
					console.log('      j:', j);

					const p = vertices[i];
					const q = vertices[j];
					console.log('      p:', p);
					console.log('      q:', q);

					const diff = math.subtract(p, q);
					console.log('      diff:', diff);

					let diff_norm = math.norm(diff);
					console.log('      diff_norm (before clamping):', diff_norm);

					// Handle near-zero distances
					if (diff_norm < 1e-8) {
						diff_norm = 1e-8;
						console.log('      diff_norm (after clamping):', diff_norm);
					}

					const term1 = 1 / Math.pow(diff_norm, 2 * sigma + 1);
					console.log('      term1:', term1);
					if (!Number.isFinite(term1)) {
						console.error(
							'      term1 is not finite! diff_norm:',
							diff_norm,
							'2 * sigma + 1:',
							2 * sigma + 1
						);
					}
					elt1 += term1;
					console.log('      elt1 (after adding term1):', elt1);

					const alph = 2; // For B0
					const bet = 4; // For B0
					console.log('      alph (for B0):', alph);
					console.log('      bet (for B0):', bet);

					// det diff, edgeTangents[I]
					const cross = math.det([diff, edgeTangents[I]]);
					console.log('      cross product:', cross);

					const cross_norm = math.norm(cross);
					console.log('      cross_norm:', cross_norm);

					const k_numerator = Math.pow(cross_norm, alph);
					console.log('      k_numerator:', k_numerator);

					const k_denominator = Math.pow(diff_norm, bet);
					console.log('      k_denominator:', k_denominator);

					const k = k_numerator / k_denominator;
					console.log('      k:', k);
					if (!Number.isFinite(k)) {
						console.error(
							'      k is not finite! k_numerator:',
							k_numerator,
							'k_denominator:',
							k_denominator
						);
					}

					const term2 = k / Math.pow(diff_norm, 2 * sigma + 1);
					console.log('      term2:', term2);
					if (!Number.isFinite(term2)) {
						console.error(
							'      term2 is not finite! k:',
							k,
							'diff_norm:',
							diff_norm,
							'2 * sigma + 1:',
							2 * sigma + 1
						);
					}
					elt2 += term2;
					console.log('      elt2 (after adding term2):', elt2);
				}
			}
			console.log('  elt1 (final):', elt1);
			console.log('  elt2 (final):', elt2);

			const w_ij_factor = 0.25 * edgeLengths[I] * edgeLengths[J];
			console.log('  w_ij_factor:', w_ij_factor);

			const W_IJ = w_ij_factor * elt1;
			console.log('  W_IJ:', W_IJ);
			if (!Number.isFinite(W_IJ)) {
				console.error('  W_IJ is not finite! w_ij_factor:', w_ij_factor, 'elt1:', elt1);
			}
			W.set([I, J], W_IJ);
			console.log('  W (after setting):', W);

			const W0_IJ = w_ij_factor * elt2;
			console.log('  W0_IJ:', W0_IJ);
			if (!Number.isFinite(W0_IJ)) {
				console.error('  W0_IJ is not finite! w_ij_factor:', w_ij_factor, 'elt2:', elt2);
			}
			W0.set([I, J], W0_IJ);
			console.log('  W0 (after setting):', W0);
		}
	}

	console.log('Returning from build_weights');
	return { W, W0 };
}

function calculateLowOrderTerm(vertices, edges, sigma, W0) {
	const numVertices = vertices.length;
	const B0 = math.zeros(numVertices, numVertices);
	const disjointEdges = calculateDisjointEdgePairs(edges); // Use your existing function

	for (let I = 0; I < edges.length; I++) {
		for (const J of disjointEdges[I]) {
			const w_IJ_0 = W0.get([I, J]); // Get the precomputed weight

			// Apply the increments
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

function calculateHighOrderTerm(vertices, edges, sigma, W) {
	const numVertices = vertices.length;
	const B = math.zeros(numVertices, numVertices);
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
	const disjointEdges = calculateDisjointEdgePairs(edges);

	for (let I = 0; I < edges.length; I++) {
		for (const J of disjointEdges[I]) {
			const l_I = edgeLengths[I];
			const l_J = edgeLengths[J];
			const T_I = edgeTangents[I];
			const T_J = edgeTangents[J];

			const w_IJ = W.get([I, J]); // Get the precomputed weight
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

export function calculateDiscreteInnerProduct(vertices, edges, edgeTangents, alpha, beta, B, B0) {
	console.log('Calculating discrete inner product with alpha:', alpha, 'beta:', beta);
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

export function build_A_bar_2D(alpha, beta, vertices, edges, sigma) {
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
	const B = calculateHighOrderTerm(vertices, edges, sigma, W);
	const B0 = calculateLowOrderTerm(vertices, edges, sigma, W0);
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
	const sigma = (beta - 1) / alpha;
	const { A_bar, B, B0 } = build_A_bar_2D(alpha, beta, vertices, edges, sigma);
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
