// src/lib/innerProduct.js
import * as math from 'mathjs';
import { calculateEdgeProperties, calculateDisjointEdgePairs } from '$lib/energyCalculations';

function calculateWeights(vertices, edges, edgeTangents, alpha, beta, sigma) {
	console.log('Calculating weights with alpha:', alpha, 'beta:', beta, 'sigma:', sigma);
	const numEdges = edges.length;
	const disjointPairs = calculateDisjointEdgePairs(edges);
	const { edgeLengths } = calculateEdgeProperties(vertices, edges);
	const weightsHigh = math.zeros(numEdges, numEdges); // w_IJ for high-order term B
	const weightsLow = math.zeros(numEdges, numEdges); // w_IJ^0 for low-order term B^0

	const s = (beta - 1) / alpha; // Paper's s = (β - 1) / α
	const sigmaCorrected = s - 1; // σ = s - 1 as per paper

	for (let I = 0; I < numEdges; I++) {
		for (const J of disjointPairs[I]) {
			const x_I = [
				(vertices[edges[I][0]][0] + vertices[edges[I][1]][0]) / 2,
				(vertices[edges[I][0]][1] + vertices[edges[I][1]][1]) / 2
			];
			const x_J = [
				(vertices[edges[J][0]][0] + vertices[edges[J][1]][0]) / 2,
				(vertices[edges[J][0]][1] + vertices[edges[J][1]][1]) / 2
			];

			const diff = math.subtract(x_I, x_J);
			const dist = math.norm(diff) + 1e-6;

			// High-order weight w_IJ (Equation 25, simplified for discrete edges)
			const kappaHigh = 2 * sigmaCorrected + 1; // For high-order term B
			weightsHigh.set(
				[I, J],
				(1 / (4 * edgeLengths[I] * edgeLengths[J])) * Math.pow(dist, -kappaHigh)
			);
			weightsHigh.set([J, I], weightsHigh.get([I, J])); // Symmetric

			// Low-order weight w_IJ^0 (Equation 33, includes k_4^2 term approximated as cross product squared)
			const kappaLow = 2 * sigmaCorrected + 1; // For low-order term B^0
			const T_I = edgeTangents[I];
			const cross_IJ = T_I[0] * diff[1] - T_I[1] * diff[0];
			const k4Squared = Math.pow(cross_IJ, 2); // Simplified k_4^2 for 2D
			weightsLow.set(
				[I, J],
				(1 / (4 * edgeLengths[I] * edgeLengths[J])) * (k4Squared / Math.pow(dist, kappaLow))
			);
			weightsLow.set([J, I], weightsLow.get([I, J])); // Symmetric
		}
	}

	console.log('High-order weights:', weightsHigh.toArray());
	console.log('Low-order weights:', weightsLow.toArray());
	return { weightsHigh, weightsLow };
}

function calculateLowOrderTerm(vertices, edges, edgeTangents, alpha, beta) {
	console.log('Computing B0 - Low-order term');
	const numEdges = edges.length;
	const numVertices = vertices.length;
	const { edgeLengths } = calculateEdgeProperties(vertices, edges);
	const s = (beta - 1) / alpha;
	const sigma = s - 1; // σ = s - 1
	const { weightsLow } = calculateWeights(vertices, edges, edgeTangents, alpha, beta, sigma);
	const B0 = math.zeros(numVertices * 2, numVertices * 2); // 2D vertex-based matrix

	for (let I = 0; I < numEdges; I++) {
		const i1 = edges[I][0]; // Start vertex of edge I
		const i2 = edges[I][1]; // End vertex of edge I
		const x_I = [(vertices[i1][0] + vertices[i2][0]) / 2, (vertices[i1][1] + vertices[i2][1]) / 2];
		const T_I = edgeTangents[I];
		console.log(`B0 Edge ${I}: x_I = ${x_I}, T_I = ${T_I}`);

		for (const J of calculateDisjointEdgePairs(edges)[I]) {
			const j1 = edges[J][0];
			const j2 = edges[J][1];
			const x_J = [
				(vertices[j1][0] + vertices[j2][0]) / 2,
				(vertices[j1][1] + vertices[j2][1]) / 2
			];

			const diff = math.subtract(x_I, x_J);
			const weight = weightsLow.get([I, J]);

			// Low-order term increments (Equation 33 and increments)
			for (let a = 0; a < 2; a++) {
				// x, y coordinates (1, 2 in paper)
				for (let b = 0; b < 2; b++) {
					// i1,i1 increment
					B0.set([i1 * 2 + a, i1 * 2 + b], B0.get([i1 * 2 + a, i1 * 2 + b]) + (1 / 4) * weight);
					// i1,j1 increment
					B0.set([i1 * 2 + a, j1 * 2 + b], B0.get([i1 * 2 + a, j1 * 2 + b]) - (1 / 4) * weight);
					// j1,i1 increment
					B0.set([j1 * 2 + a, i1 * 2 + b], B0.get([j1 * 2 + a, i1 * 2 + b]) - (1 / 4) * weight);
					// j1,j1 increment
					B0.set([j1 * 2 + a, j1 * 2 + b], B0.get([j1 * 2 + a, j1 * 2 + b]) + (1 / 4) * weight);
				}
			}
		}
	}

	console.log('B0 Matrix (vertex-based):', B0.toArray());
	return B0;
}

function calculateHighOrderTerm(vertices, edges, edgeTangents, alpha, beta) {
	console.log('Computing B - High-order term');
	const numEdges = edges.length;
	const numVertices = vertices.length;
	const { edgeLengths } = calculateEdgeProperties(vertices, edges);
	const s = (beta - 1) / alpha;
	const sigma = s - 1; // σ = s - 1
	const { weightsHigh } = calculateWeights(vertices, edges, edgeTangents, alpha, beta, sigma);
	const B = math.zeros(numVertices * 2, numVertices * 2); // 2D vertex-based matrix

	for (let I = 0; I < numEdges; I++) {
		const i1 = edges[I][0]; // Start vertex
		const i2 = edges[I][1]; // End vertex
		const x_I = [(vertices[i1][0] + vertices[i2][0]) / 2, (vertices[i1][1] + vertices[i2][1]) / 2];
		const T_I = edgeTangents[I];
		console.log(`B Edge ${I}: x_I = ${x_I}, T_I = ${T_I}`);

		for (const J of calculateDisjointEdgePairs(edges)[I]) {
			const j1 = edges[J][0];
			const j2 = edges[J][1];
			const x_J = [
				(vertices[j1][0] + vertices[j2][0]) / 2,
				(vertices[j1][1] + vertices[j2][1]) / 2
			];
			const T_J = edgeTangents[J];

			const diff = math.subtract(x_I, x_J);
			const dist = math.norm(diff) + 1e-6;
			const weight = weightsHigh.get([I, J]);
			const dotProduct = T_I[0] * T_J[0] + T_I[1] * T_J[1]; // <T_I, T_J>

			// High-order term increments (Equation 21 and increments)
			for (let a = 0; a < 2; a++) {
				// x, y coordinates (1, 2 in paper)
				for (let b = 0; b < 2; b++) {
					const sign = (-1) ** (a + b);
					// i1,i1 increment
					B.set(
						[i1 * 2 + a, i1 * 2 + b],
						B.get([i1 * 2 + a, i1 * 2 + b]) + (sign * weight) / (edgeLengths[I] * edgeLengths[I])
					);
					// i1,j1 increment
					B.set(
						[i1 * 2 + a, j1 * 2 + b],
						B.get([i1 * 2 + a, j1 * 2 + b]) -
							(sign * weight * dotProduct) / (edgeLengths[I] * edgeLengths[J])
					);
					// j1,j1 increment
					B.set(
						[j1 * 2 + a, j1 * 2 + b],
						B.get([j1 * 2 + a, j1 * 2 + b]) + (sign * weight) / (edgeLengths[J] * edgeLengths[J])
					);
					// j1,i1 increment
					B.set(
						[j1 * 2 + a, i1 * 2 + b],
						B.get([j1 * 2 + a, i1 * 2 + b]) -
							(sign * weight * dotProduct) / (edgeLengths[J] * edgeLengths[I])
					);
				}
			}
		}
	}

	console.log('B Matrix (vertex-based):', B.toArray());
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
	console.log('Computing preconditioned gradient with alpha:', alpha, 'beta:', beta);
	console.log('Input vertices:', vertices);
	console.log('Input L2Gradient:', L2Gradient);
	const { A } = calculateDiscreteInnerProduct(vertices, edges, edgeTangents, alpha, beta);
	const numVertices = vertices.length;
	const dim = 2;

	// Convert L2Gradient from vertex-based to flat form
	const gradientFlat = [];
	for (let i = 0; i < numVertices; i++) {
		gradientFlat.push(L2Gradient[i][0], L2Gradient[i][1]);
	}
	console.log('GradientFlat:', gradientFlat);

	const AFull = A; // Already vertex-based and 2D
	// Add regularization to ensure invertibility
	const epsilon = 1e-2;
	for (let i = 0; i < numVertices * dim; i++) {
		AFull.set([i, i], AFull.get([i, i]) + epsilon);
	}
	console.log('AFull:', AFull.toArray());
	console.log('AFull determinant:', math.det(AFull));

	const negGradientFlat = math.multiply(-1, gradientFlat);
	console.log('NegGradientFlat:', negGradientFlat);

	let gradFlatFull;
	try {
		gradFlatFull = math.multiply(math.inv(AFull), negGradientFlat);
		console.log('gradFlatFull after solve:', gradFlatFull.toArray());
	} catch (e) {
		console.error('Matrix inversion failed:', e);
		gradFlatFull = negGradientFlat; // Fallback to L2 gradient
	}

	const gradFlatArray = Array.isArray(gradFlatFull) ? gradFlatFull : gradFlatFull.toArray();
	const preconditionedGradient = [];
	for (let i = 0; i < numVertices; i++) {
		const x = isNaN(gradFlatArray[i * dim]) ? 0 : gradFlatArray[i * dim];
		const y = isNaN(gradFlatArray[i * dim + 1]) ? 0 : gradFlatArray[i * dim + 1];
		preconditionedGradient[i] = [x, y];
		console.log(`PreconditionedGradient[${i}]: [${x}, ${y}]`);
	}
	return preconditionedGradient;
}
