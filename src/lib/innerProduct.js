// src/lib/innerProduct.js
import * as math from 'mathjs';
import {
	calculateEdgeProperties,
	calculateDisjointEdgePairs,
	tangentPointKernel
} from '$lib/energyCalculations';

function calculateLowOrderTerm(vertices, edges, sigma, alpha = 2, beta = 4) {
	alpha = 2;
	beta = 4; // hardcode it for now, as per paper.
	const numVertices = vertices.length;
	const B0 = math.zeros(numVertices, numVertices);
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
	const disjointEdges = calculateDisjointEdgePairs(edges); // Use your existing function

	for (let I = 0; I < edges.length; I++) {
		for (const J of disjointEdges[I]) {
			// Iterate directly over disjoint edges of I
			const [i1, i2] = edges[I];
			const [j1, j2] = edges[J];
			const l_I = edgeLengths[I];
			const l_J = edgeLengths[J];
			const T_I = edgeTangents[I];

			let w_IJ_0 = 0;
			const pairs = [
				[i1, j1],
				[i1, j2],
				[i2, j1],
				[i2, j2]
			];
			for (const [i, j] of pairs) {
				const diff = math.subtract(vertices[i], vertices[j]);
				const dist = math.norm(diff) + 1e-8;
				const kernel = tangentPointKernel(vertices[i], vertices[j], T_I, alpha, beta); // Use your kernel function
				w_IJ_0 += kernel / Math.pow(dist, 2 * sigma + 1 - beta); // Correct denominator
			}
			w_IJ_0 *= (l_I * l_J) / 4;

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

	return B0;
}

function calculateHighOrderTerm(vertices, edges, sigma) {
	const numVertices = vertices.length;
	const B = math.zeros(numVertices, numVertices);
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
	const disjointEdges = calculateDisjointEdgePairs(edges); // Use your existing function

	for (let I = 0; I < edges.length; I++) {
		for (const J of disjointEdges[I]) {
			// Iterate directly over disjoint edges of I
			const [i1, i2] = edges[I];
			const [j1, j2] = edges[J];
			const l_I = edgeLengths[I];
			const l_J = edgeLengths[J];
			const T_I = edgeTangents[I];
			const T_J = edgeTangents[J];

			let w_IJ = 0;
			const pairs = [
				[i1, j1],
				[i1, j2],
				[i2, j1],
				[i2, j2]
			];
			for (const [i, j] of pairs) {
				const diff = math.subtract(vertices[i], vertices[j]);
				const dist = math.norm(diff) + 1e-8; // Add small epsilon for numerical stability
				w_IJ += 1 / Math.pow(dist, 2 * sigma + 1);
			}
			w_IJ *= (l_I * l_J) / 4;

			const dot_TI_TJ = math.dot(T_I, T_J);

			// Apply the increments as described in the paper
			for (let a = 0; a < 2; a++) {
				for (let b = 0; b < 2; b++) {
					const sign = Math.pow(-1, a + b);
					const i_a = edges[I][a];
					const i_b = edges[I][b];
					const j_a = edges[J][a];
					const j_b = edges[J][b];

					B.set([i_a, i_b], B.get([i_a, i_b]) + (sign * w_IJ) / (l_I * l_I));
					B.set([j_a, j_b], B.get([j_a, j_b]) + (sign * w_IJ) / (l_J * l_J));
					B.set([i_a, j_b], B.get([i_a, j_b]) - (sign * w_IJ * dot_TI_TJ) / (l_I * l_J));
					B.set([j_a, i_b], B.get([j_a, i_b]) - (sign * w_IJ * dot_TI_TJ) / (l_I * l_J));
				}
			}
		}
	}

	return B;
}

export function calculateDiscreteInnerProduct(vertices, edges, edgeTangents, alpha, beta) {
	console.log('Calculating discrete inner product with alpha:', alpha, 'beta:', beta);
	const s = (beta - 1) / alpha;
	const sigma = s - 1;
	console.log('Computed s:', s, 'Sigma (s - 1):', sigma);

	const B0 = calculateLowOrderTerm(vertices, edges, edgeTangents, sigma, alpha, beta); // Pass alpha and beta
	const B = calculateHighOrderTerm(vertices, edges, edgeTangents, sigma);
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

export function build_A_bar_2D(vertices, edges, sigma) {
	const B = calculateHighOrderTerm(vertices, edges, sigma);
	const B0 = calculateLowOrderTerm(vertices, edges, sigma);
	const A = math.add(B, B0);

	const numVertices = vertices.length;
	const A_bar = math.zeros(2 * numVertices, 2 * numVertices); // 2D, so 2 * numVertices

	// Block-diagonal construction (only two blocks for 2D)
	A_bar.subset(math.index(math.range(0, numVertices), math.range(0, numVertices)), A);
	A_bar.subset(
		math.index(math.range(numVertices, 2 * numVertices), math.range(numVertices, 2 * numVertices)),
		A
	);

	return A_bar;
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
	const { A } = calculateDiscreteInnerProduct(vertices, edges, edgeTangents, alpha, beta);
	const numVertices = vertices.length;
	const AFull = math.zeros(numVertices * 2, numVertices * 2);

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
