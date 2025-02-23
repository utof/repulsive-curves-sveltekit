// src/lib/innerProduct.js
import * as math from 'mathjs';
import { calculateEdgeProperties, calculateDisjointEdgePairs } from '$lib/energyCalculations';

// function calculateWeights(vertices, edges, edgeTangents, alpha, beta, sigma) {
// 	console.log('Calculating weights with alpha:', alpha, 'beta:', beta, 'sigma:', sigma);
// 	const numEdges = edges.length;
// 	const disjointPairs = calculateDisjointEdgePairs(edges);
// 	const { edgeLengths } = calculateEdgeProperties(vertices, edges);
// 	const weightsHigh = math.zeros(numEdges, numEdges); // w_IJ for high-order term B
// 	const weightsLow = math.zeros(numEdges, numEdges); // w_IJ^0 for low-order term B^0

// 	const s = (beta - 1) / alpha; // Paper's s = (β - 1) / α
// 	const sigmaCorrected = s - 1; // σ = s - 1 as per paper

// 	for (let I = 0; I < numEdges; I++) {
// 		for (const J of disjointPairs[I]) {
// 			const x_I = [
// 				(vertices[edges[I][0]][0] + vertices[edges[I][1]][0]) / 2,
// 				(vertices[edges[I][0]][1] + vertices[edges[I][1]][1]) / 2
// 			];
// 			const x_J = [
// 				(vertices[edges[J][0]][0] + vertices[edges[J][1]][0]) / 2,
// 				(vertices[edges[J][0]][1] + vertices[edges[J][1]][1]) / 2
// 			];

// 			const diff = math.subtract(x_I, x_J);
// 			const dist = math.norm(diff) + 1e-6;

// 			// High-order weight w_IJ (Equation 25, simplified for discrete edges)
// 			const kappaHigh = 2 * sigmaCorrected + 1; // For high-order term B
// 			weightsHigh.set(
// 				[I, J],
// 				(1 / (4 * edgeLengths[I] * edgeLengths[J])) * Math.pow(dist, -kappaHigh)
// 			);
// 			weightsHigh.set([J, I], weightsHigh.get([I, J])); // Symmetric

// 			// Low-order weight w_IJ^0 (Equation 33, includes k_4^2 term approximated as cross product squared)
// 			const kappaLow = 2 * sigmaCorrected + 1; // For low-order term B^0
// 			const T_I = edgeTangents[I];
// 			const cross_IJ = T_I[0] * diff[1] - T_I[1] * diff[0];
// 			const k4Squared = Math.pow(cross_IJ, 2); // Simplified k_4^2 for 2D
// 			weightsLow.set(
// 				[I, J],
// 				(1 / (4 * edgeLengths[I] * edgeLengths[J])) * (k4Squared / Math.pow(dist, kappaLow))
// 			);
// 			weightsLow.set([J, I], weightsLow.get([I, J])); // Symmetric
// 		}
// 	}

// 	console.log('High-order weights:', weightsHigh.toArray());
// 	console.log('Low-order weights:', weightsLow.toArray());
// 	return { weightsHigh, weightsLow };
// }

function calculateLowOrderTerm(vertices, edges, edgeTangents, sigma) {
	const { edgeLengths } = calculateEdgeProperties(vertices, edges);
	const numVertices = vertices.length;
	const B0 = math.zeros(numVertices, numVertices);
	const kappa = 2 * sigma + 1; // Match energy’s exponent - 4
	for (let I = 0; I < edges.length; I++) {
		const i1 = edges[I][0],
			i2 = edges[I][1];
		const l_I = edgeLengths[I];
		const T_I = edgeTangents[I];
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
			for (const [i, j] of pairs) {
				const diff = math.subtract(vertices[i], vertices[j]);
				const dist = math.norm(diff) + 1e-6;
				const cross = T_I[0] * diff[1] - T_I[1] * diff[0];
				w_IJ += Math.pow(Math.abs(cross), 2) / Math.pow(dist, kappa + 4);
			}
			w_IJ *= (l_I * l_J) / 4;
			B0.set([i1, i1], B0.get([i1, i1]) + w_IJ);
			B0.set([i2, i2], B0.get([i2, i2]) + w_IJ);
			B0.set([i1, i2], B0.get([i1, i2]) - w_IJ);
			B0.set([i2, i1], B0.get([i2, i1]) - w_IJ);
		}
	}
	return B0;
}

function calculateHighOrderTerm(vertices, edges, edgeTangents, sigma) {
	const { edgeLengths } = calculateEdgeProperties(vertices, edges);
	const numVertices = vertices.length;
	const B = math.zeros(numVertices, numVertices);
	const kappa = 2 * sigma + 1;
	for (let I = 0; I < edges.length; I++) {
		const i1 = edges[I][0],
			i2 = edges[I][1];
		const l_I = edgeLengths[I];
		const T_I = edgeTangents[I];
		for (let J = 0; J < edges.length; J++) {
			if (I === J || edges[I][0] === edges[J][0] || edges[I][1] === edges[J][1]) continue;
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
			for (const [i, j] of pairs) {
				const diff = math.subtract(vertices[i], vertices[j]);
				const dist = math.norm(diff) + 1e-6;
				w_IJ += 1 / Math.pow(dist, kappa);
			}
			w_IJ *= (l_I * l_J) / 4;
			const dotTT = math.dot(T_I, T_J);
			B.set([i1, i1], B.get([i1, i1]) + w_IJ / (l_I * l_I));
			B.set([i2, i2], B.get([i2, i2]) + w_IJ / (l_I * l_I));
			B.set([i1, i2], B.get([i1, i2]) - w_IJ / (l_I * l_I));
			B.set([i2, i1], B.get([i2, i1]) - w_IJ / (l_I * l_I));
			// Cross terms with J
			B.set([j1, j1], B.get([j1, j1]) + w_IJ / (l_J * l_J));
			B.set([j2, j2], B.get([j2, j2]) + w_IJ / (l_J * l_J));
			B.set([j1, j2], B.get([j1, j2]) - w_IJ / (l_J * l_J));
			B.set([j2, j1], B.get([j2, j1]) - w_IJ / (l_J * l_J));
			// Interaction terms
			B.set([i1, j1], B.get([i1, j1]) - (w_IJ * dotTT) / (l_I * l_J));
			B.set([i1, j2], B.get([i1, j2]) + (w_IJ * dotTT) / (l_I * l_J));
			B.set([i2, j1], B.get([i2, j1]) + (w_IJ * dotTT) / (l_I * l_J));
			B.set([i2, j2], B.get([i2, j2]) - (w_IJ * dotTT) / (l_I * l_J));
		}
	}
	return B;
}

export function calculateDiscreteInnerProduct(vertices, edges, edgeTangents, alpha, beta) {
	console.log('Calculating discrete inner product with alpha:', alpha, 'beta:', beta);
	const s = (beta - 1) / alpha;
	const sigma = s - 1; // σ = s - 1
	console.log('Sigma (s - 1):', sigma);
	const B0 = calculateLowOrderTerm(vertices, edges, edgeTangents, alpha, beta);
	const B = calculateHighOrderTerm(vertices, edges, edgeTangents, alpha, beta);
	const A = math.add(B0, B);
	console.log('A Matrix:', A.toArray());
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
	const { A, B0, B } = calculateDiscreteInnerProduct(vertices, edges, edgeTangents, alpha, beta);
	const numVertices = vertices.length;
	const AFull = math.zeros(numVertices * 2, numVertices * 2);
	for (let i = 0; i < numVertices; i++) {
		for (let j = 0; j < numVertices; j++) {
			const val = A.get([i, j]);
			AFull.set([i * 2, j * 2], val);
			AFull.set([i * 2 + 1, j * 2 + 1], val);
		}
	}
	const gradientFlat = L2Gradient.flat();
	const negGradientFlat = math.multiply(-1, gradientFlat);
	let gradFlatFull;
	try {
		gradFlatFull = math.multiply(math.inv(AFull), negGradientFlat);
	} catch (e) {
		console.error('Inversion failed:', e);
		gradFlatFull = math.zeros(numVertices * 2);
	}
	const grad = [];
	for (let i = 0; i < numVertices; i++) {
		grad[i] = [gradFlatFull[i * 2], gradFlatFull[i * 2 + 1]];
	}
	return grad;
}
