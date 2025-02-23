// src/lib/innerProduct.js
import * as math from 'mathjs';
import { calculateEdgeProperties } from '$lib/energyCalculations';

function calculateWeights(edgeLengths) {
	console.log('Calculating weights:', edgeLengths);
	return edgeLengths.map((length) => {
		const weight = isNaN(length) ? 0 : length;
		console.log(`Weight for length ${length}: ${weight}`);
		return weight;
	});
}

function calculateLowOrderTerm(vertices, edges, edgeTangents, sigma) {
	console.log('Computing B0 - Low-order term');
	const { edgeLengths } = calculateEdgeProperties(vertices, edges);
	const numEdges = edges.length;
	const weights = calculateWeights(edgeLengths);
	const B0 = math.zeros(numEdges, numEdges);

	const kappa = 2 * sigma + 5;
	for (let I = 0; I < numEdges; I++) {
		const x_I = [
			(vertices[edges[I][0]][0] + vertices[edges[I][1]][0]) / 2,
			(vertices[edges[I][0]][1] + vertices[edges[I][1]][1]) / 2
		];
		const T_I = edgeTangents[I];
		console.log(`B0 Edge ${I}: x_I = ${x_I}, T_I = ${T_I}`);

		for (let J = 0; J < numEdges; J++) {
			if (I === J) continue;

			const x_J = [
				(vertices[edges[J][0]][0] + vertices[edges[J][1]][0]) / 2,
				(vertices[edges[J][0]][1] + vertices[edges[J][1]][1]) / 2
			];
			const T_J = edgeTangents[J];

			const diff = math.subtract(x_I, x_J);
			const dist = math.norm(diff) + 1e-6;
			const cross_I = T_I[0] * diff[1] - T_I[1] * diff[0];
			const cross_J = T_J[0] * diff[1] - T_J[1] * diff[0];

			const term_I = Math.pow(Math.abs(cross_I), 2) / Math.pow(dist, kappa);
			const term_J = Math.pow(Math.abs(cross_J), 2) / Math.pow(dist, kappa);
			const value = (term_I + term_J) * weights[I] * weights[J];

			console.log(`B0 [${I},${J}]: diff = ${diff}, dist = ${dist}, value = ${value}`);
			B0.set([I, J], -value);
		}

		let diagSum = 0;
		for (let J = 0; J < numEdges; J++) {
			if (I !== J) diagSum += Math.abs(B0.get([I, J]));
		}
		B0.set([I, I], diagSum);
		console.log(`B0 Diagonal [${I},${I}]: ${diagSum}`);
	}
	console.log('B0 Matrix:', B0.toArray());
	return B0;
}

function calculateHighOrderTerm(vertices, edges, edgeTangents, sigma) {
	console.log('Computing B - High-order term');
	const { edgeLengths } = calculateEdgeProperties(vertices, edges);
	const numEdges = edges.length;
	const weights = calculateWeights(edgeLengths);
	const B = math.zeros(numEdges, numEdges);

	const kappa = 2 * sigma + 1;
	for (let I = 0; I < numEdges; I++) {
		const x_I = [
			(vertices[edges[I][0]][0] + vertices[edges[I][1]][0]) / 2,
			(vertices[edges[I][0]][1] + vertices[edges[I][1]][1]) / 2
		];
		const T_I = edgeTangents[I];
		console.log(`B Edge ${I}: x_I = ${x_I}, T_I = ${T_I}`);

		for (let J = 0; J < numEdges; J++) {
			if (I === J) continue;

			const x_J = [
				(vertices[edges[J][0]][0] + vertices[edges[J][1]][0]) / 2,
				(vertices[edges[J][0]][1] + vertices[edges[J][1]][1]) / 2
			];
			const T_J = edgeTangents[J];

			const diff = math.subtract(x_I, x_J);
			const dist = math.norm(diff) + 1e-6;
			const cross_I = T_I[0] * diff[1] - T_I[1] * diff[0];
			const cross_J = T_J[0] * diff[1] - T_J[1] * diff[0];

			const term_I = Math.pow(Math.abs(cross_I), 2) / Math.pow(dist, kappa);
			const term_J = Math.pow(Math.abs(cross_J), 2) / Math.pow(dist, kappa);
			const value = (term_I + term_J) * weights[I] * weights[J];

			console.log(`B [${I},${J}]: diff = ${diff}, dist = ${dist}, value = ${value}`);
			B.set([I, J], -value);
		}

		let diagSum = 0;
		for (let J = 0; J < numEdges; J++) {
			if (I !== J) diagSum += Math.abs(B.get([I, J]));
		}
		B.set([I, I], diagSum);
		console.log(`B Diagonal [${I},${I}]: ${diagSum}`);
	}
	console.log('B Matrix:', B.toArray());
	return B;
}

export function calculateDiscreteInnerProduct(vertices, edges, edgeTangents, sigma) {
	console.log('Calculating discrete inner product');
	const B0 = calculateLowOrderTerm(vertices, edges, edgeTangents, sigma);
	const B = calculateHighOrderTerm(vertices, edges, edgeTangents, sigma);
	const A = math.add(B0, B);
	console.log('A Matrix:', A.toArray());
	return { A, B0, B };
}

export function computePreconditionedGradient(vertices, edges, edgeTangents, sigma, L2Gradient) {
	console.log('Computing preconditioned gradient');
	console.log('Input vertices:', vertices);
	console.log('Input L2Gradient:', L2Gradient);
	const { A } = calculateDiscreteInnerProduct(vertices, edges, edgeTangents, sigma);
	const numEdges = edges.length;

	const gradientFlat = [];
	for (let I = 0; I < numEdges; I++) {
		const v0 = edges[I][0];
		const v1 = edges[I][1];
		const avgX = (L2Gradient[v0][0] + L2Gradient[v1][0]) / 2;
		const avgY = (L2Gradient[v0][1] + L2Gradient[v1][1]) / 2;
		gradientFlat.push(isNaN(avgX) ? 0 : avgX, isNaN(avgY) ? 0 : avgY);
		console.log(`GradientFlat[${I}]: [${avgX}, ${avgY}]`);
	}

	const dim = 2;
	const AFull = math.zeros(numEdges * dim, numEdges * dim);
	for (let i = 0; i < numEdges; i++) {
		for (let j = 0; j < numEdges; j++) {
			const val = A.get([i, j]);
			AFull.set([i * dim, j * dim], val);
			AFull.set([i * dim + 1, j * dim + 1], val);
		}
	}
	// Add regularization to ensure invertibility
	const epsilon = 1e-6;
	for (let i = 0; i < numEdges * dim; i++) {
		AFull.set([i, i], AFull.get([i, i]) + epsilon);
	}
	console.log('AFull Matrix (regularized):', AFull.toArray());
	console.log('AFull determinant:', math.det(AFull));

	const negGradientFlat = math.multiply(-1, gradientFlat);
	console.log('NegGradientFlat:', negGradientFlat);

	let gradFlatFull;
	try {
		gradFlatFull = math.multiply(math.inv(AFull), negGradientFlat);
		console.log('gradFlatFull after solve (raw):', gradFlatFull.toArray());
	} catch (e) {
		console.error('Matrix inversion failed:', e);
		gradFlatFull = math.zeros(numEdges * dim); // Fallback
	}

	const gradFlatArray = Array.isArray(gradFlatFull) ? gradFlatFull : gradFlatFull.toArray();

	const gradFlat = [];
	for (let i = 0; i < numEdges; i++) {
		const x = isNaN(gradFlatArray[i * dim]) ? 0 : gradFlatArray[i * dim];
		const y = isNaN(gradFlatArray[i * dim + 1]) ? 0 : gradFlatArray[i * dim + 1];
		gradFlat[i] = [x, y];
		console.log(`gradFlat[${i}]: [${x}, ${y}]`);
	}

	const preconditionedGradient = vertices.map(() => [0, 0]);
	const vertexCounts = vertices.map(() => 0);
	for (let I = 0; I < numEdges; I++) {
		const v0 = edges[I][0];
		const v1 = edges[I][1];
		preconditionedGradient[v0][0] += gradFlat[I][0];
		preconditionedGradient[v0][1] += gradFlat[I][1];
		preconditionedGradient[v1][0] += gradFlat[I][0];
		preconditionedGradient[v1][1] += gradFlat[I][1];
		vertexCounts[v0]++;
		vertexCounts[v1]++;
	}

	const result = preconditionedGradient.map((grad, i) => {
		const x = grad[0] / (vertexCounts[i] || 1);
		const y = grad[1] / (vertexCounts[i] || 1);
		console.log(`PreconditionedGradient[${i}]: [${x}, ${y}]`);
		return [x, y];
	});
	return result;
}
