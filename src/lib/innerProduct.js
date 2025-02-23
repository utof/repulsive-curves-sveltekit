// src/lib/innerProduct.js
import * as math from 'mathjs';
import { calculateEdgeProperties, calculateDisjointEdgePairs } from '$lib/energyCalculations';

function calculateLowOrderTerm(vertices, edges, edgeTangents, sigma) {
	console.log('Starting calculateLowOrderTerm with sigma:', sigma);
	const { edgeLengths } = calculateEdgeProperties(vertices, edges);
	const numVertices = vertices.length;
	const B0 = math.zeros(numVertices, numVertices);
	const kappa = 2 * sigma + 1;
	console.log('kappa for B0:', kappa);

	for (let I = 0; I < edges.length; I++) {
		const i1 = edges[I][0],
			i2 = edges[I][1];
		const l_I = edgeLengths[I];
		const T_I = edgeTangents[I];
		console.log(`B0 Edge ${I}: [${i1}, ${i2}], length: ${l_I}, tangent: [${T_I}]`);

		for (let J = 0; J < edges.length; J++) {
			if (I === J || edges[I][0] === edges[J][0] || edges[I][1] === edges[J][1]) continue;
			const j1 = edges[J][0],
				j2 = edges[J][1];
			const l_J = edgeLengths[J];
			let w_IJ = 0;
			const pairs = [
				[i1, j1],
				[i1, j2],
				[i2, j1],
				[i2, j2]
			];
			console.log(`  Pair with Edge ${J}: [${j1}, ${j2}]`);

			for (const [i, j] of pairs) {
				const diff = math.subtract(vertices[i], vertices[j]);
				const dist = math.norm(diff) + 1e-6;
				const cross = T_I[0] * diff[1] - T_I[1] * diff[0];
				const contrib = Math.pow(Math.abs(cross), 2) / Math.pow(dist, kappa);
				w_IJ += contrib;
				console.log(`    Pair [${i}, ${j}]: dist=${dist}, cross=${cross}, contrib=${contrib}`);
			}
			w_IJ *= (l_I * l_J) / 4;
			console.log(`  w_IJ for ${I}-${J}: ${w_IJ}`);

			B0.set([i1, i1], B0.get([i1, i1]) + w_IJ);
			B0.set([i2, i2], B0.get([i2, i2]) + w_IJ);
			B0.set([i1, i2], B0.get([i1, i2]) - w_IJ);
			B0.set([i2, i1], B0.get([i2, i1]) - w_IJ);
		}
	}
	console.log('B0 Matrix:', B0.toArray());
	return B0;
}

function calculateHighOrderTerm(vertices, edges, edgeTangents, sigma) {
	console.log('Starting calculateHighOrderTerm with sigma:', sigma);
	const { edgeLengths } = calculateEdgeProperties(vertices, edges);
	const numVertices = vertices.length;
	const B = math.zeros(numVertices, numVertices);
	const kappa = 2 * sigma + 1;
	console.log('kappa for B:', kappa);

	for (let I = 0; I < edges.length; I++) {
		const i1 = edges[I][0],
			i2 = edges[I][1];
		const l_I = edgeLengths[I];
		const T_I = edgeTangents[I];
		console.log(`B Edge ${I}: [${i1}, ${i2}], length: ${l_I}, tangent: [${T_I}]`);

		for (let J = 0; J < edges.length; J++) {
			// Skip if edges are the same or share vertices
			if (
				I === J ||
				edges[I][0] === edges[J][0] ||
				edges[I][0] === edges[J][1] ||
				edges[I][1] === edges[J][0] ||
				edges[I][1] === edges[J][1]
			) {
				console.log(`  Skipping Edge ${J}: [${edges[J][0]}, ${edges[J][1]}] (not disjoint)`);
				continue;
			}
			const j1 = edges[J][0],
				j2 = edges[J][1];
			const l_J = edgeLengths[J];
			const T_J = edgeTangents[J];
			let w_IJ = 0;
			const pairs = [
				[i1, j1],
				[i1, j2],
				[i2, j1],
				[i2, j2]
			];
			console.log(`  Pair with Edge ${J}: [${j1}, ${j2}]`);

			for (const [i, j] of pairs) {
				const diff = math.subtract(vertices[i], vertices[j]);
				const dist = math.norm(diff) + 1e-6;
				w_IJ += 1 / Math.pow(dist, kappa);
				console.log(
					`    Pair [${i}, ${j}]: dist=${dist}, w_IJ contrib=${1 / Math.pow(dist, kappa)}`
				);
			}
			w_IJ *= (l_I * l_J) / 4;
			const dotTT = math.dot(T_I, T_J);
			console.log(`  w_IJ for ${I}-${J}: ${w_IJ}, dot(T_I, T_J): ${dotTT}`);

			B.set([i1, i1], B.get([i1, i1]) + w_IJ / (l_I * l_I));
			B.set([i2, i2], B.get([i2, i2]) + w_IJ / (l_I * l_I));
			B.set([i1, i2], B.get([i1, i2]) - w_IJ / (l_I * l_I));
			B.set([i2, i1], B.get([i2, i1]) - w_IJ / (l_I * l_I));
			B.set([j1, j1], B.get([j1, j1]) + w_IJ / (l_J * l_J));
			B.set([j2, j2], B.get([j2, j2]) + w_IJ / (l_J * l_J));
			B.set([j1, j2], B.get([j1, j2]) - w_IJ / (l_J * l_J));
			B.set([j2, j1], B.get([j2, j1]) - w_IJ / (l_J * l_J));
			B.set([i1, j1], B.get([i1, j1]) - (w_IJ * dotTT) / (l_I * l_J));
			B.set([i1, j2], B.get([i1, j2]) + (w_IJ * dotTT) / (l_I * l_J));
			B.set([i2, j1], B.get([i2, j1]) + (w_IJ * dotTT) / (l_I * l_J));
			B.set([i2, j2], B.get([i2, j2]) - (w_IJ * dotTT) / (l_I * l_J));
		}
	}
	console.log('B Matrix:', B.toArray());
	return B;
}

export function calculateDiscreteInnerProduct(vertices, edges, edgeTangents, alpha, beta) {
	console.log('Calculating discrete inner product with alpha:', alpha, 'beta:', beta);
	const s = (beta - 1) / alpha;
	const sigma = s - 1;
	console.log('Computed s:', s, 'Sigma (s - 1):', sigma);

	const B0 = calculateLowOrderTerm(vertices, edges, edgeTangents, sigma);
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
	A = math.add(A, reg);
	console.log('A Matrix with regularization:', A.toArray());

	try {
		const eigenvalues = math.eigs(A).values.toArray();
		console.log('Eigenvalues of A:', eigenvalues);
		const minEigen = Math.min(...eigenvalues);
		console.log('Minimum eigenvalue:', minEigen);
		if (minEigen <= 0) console.warn('A is not positive definite!');
	} catch (e) {
		console.error('Eigenvalue computation failed:', e);
	}

	return { A, B0, B };
}

export function computePreconditionedGradient(
	vertices,
	edges,
	edgeTangents,
	alpha,
	beta,
	L2Gradient
) {
	console.log('Computing preconditioned gradient with L2Gradient:', L2Gradient);
	const { A, B0, B } = calculateDiscreteInnerProduct(vertices, edges, edgeTangents, alpha, beta);
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

	const gradientFlat = L2Gradient.flat();
	const negGradientFlat = math.multiply(-1, gradientFlat);
	console.log('Flattened negative L2 gradient:', negGradientFlat);

	let gradFlatFull;
	try {
		gradFlatFull = math.multiply(math.inv(AFull), negGradientFlat);
		console.log('Preconditioned gradient (flat) as matrix:', gradFlatFull.toArray());
	} catch (e) {
		console.error('Inversion failed:', e);
		console.log('Falling back to L2 gradient');
		gradFlatFull = math.matrix(negGradientFlat);
	}

	const gradArray = gradFlatFull.toArray();
	const grad = [];
	for (let i = 0; i < numVertices; i++) {
		grad[i] = [gradArray[i * 2], gradArray[i * 2 + 1]];
		console.log(`Gradient for vertex ${i}: [${grad[i][0]}, ${grad[i][1]}]`);
	}
	return grad;
}
