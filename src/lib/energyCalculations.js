// src/lib/energyCalculations.js
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from '$lib/stores';

let logging = false;

export function calculateEdgeProperties(vertices, edges) {
	const edgeLengths = [];
	const edgeTangents = [];
	const edgeMidpoints = [];

	for (const edge of edges) {
		const v1 = vertices[edge[0]];
		const v2 = vertices[edge[1]];

		const dx = v2[0] - v1[0];
		const dy = v2[1] - v1[1];
		const length = Math.sqrt(dx * dx + dy * dy);
		edgeLengths.push(length);

		const unitTangent = length > 0 ? [dx / length, dy / length] : [0, 0];
		edgeTangents.push(unitTangent);

		const midpoint = [
			isNaN(v1[0]) || isNaN(v2[0]) ? 0 : (v1[0] + v2[0]) / 2,
			isNaN(v1[1]) || isNaN(v2[1]) ? 0 : (v1[1] + v2[1]) / 2
		];
		edgeMidpoints.push(midpoint);

		if (logging) {
			console.log(
				`Edge [${edge[0]}, ${edge[1]}]: length = ${length}, tangent = ${unitTangent}, midpoint = ${midpoint}`
			);
		}
	}

	if (logging) {
		console.log('Edge lengths:', edgeLengths);
		console.log('Unit tangents:', edgeTangents);
		console.log('Midpoints:', edgeMidpoints);
	}

	return { edgeLengths, edgeTangents, edgeMidpoints };
}

export function tangentPointKernel(p, q, T, alpha, beta) {
	const p_ = math.matrix(p);
	const q_ = math.matrix(q);
	const T_ = math.matrix(T);
	const epsilon = get(config).epsilonKernel; // Use config epsilon

	const diff = math.subtract(p_, q_);
	const diffNorm = math.norm(diff) + epsilon; // Prevent division by zero
	const cross2D = T_.get([0]) * diff.get([1]) - T_.get([1]) * diff.get([0]); // 2D cross product (determinant)

	const numerator = Math.pow(Math.abs(cross2D), alpha);
	const denominator = Math.pow(diffNorm, beta);
	const result = numerator / denominator;

	if (logging) {
		console.log('Kernel calc:', {
			p: p_.toArray(),
			q: q_.toArray(),
			T: T_.toArray(),
			cross2D,
			diffNorm,
			numerator,
			denominator,
			result
		});
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
	console.log('Calculated disjointPairs:', disjointPairs);
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
		if (!disjointPairs[i]) {
			console.warn(`No disjoint pairs for edge ${i}`);
			continue;
		}

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
				kernelMatrix.set([j, i], sum / 4); // Keep symmetry!
			} else {
				console.warn(
					'Invalid edge index:',
					i,
					j,
					'disjointPairs:',
					disjointPairs,
					'edges.length',
					edges.length
				);
			}
		}
	}
	return kernelMatrix;
}

export function calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs) {
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
	const kernelMatrix = calculateDiscreteKernel(
		vertices,
		edges,
		edgeTangents,
		alpha,
		beta,
		disjointPairs
	);

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
	return totalEnergy / 2; // Divide by 2 because of symmetry
}

export function calculateDifferential(vertices, edges, alpha, beta, disjointPairs) {
	const h = get(config).finiteDiffH; // Use config finiteDiffH
	const numVertices = vertices.length;
	const differential = [];

	const E_original = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

	for (let i = 0; i < numVertices; i++) {
		differential[i] = [0, 0];
		for (let dim = 0; dim < 2; dim++) {
			const vertices_perturbed = vertices.map((v) => [...v]);
			vertices_perturbed[i][dim] += h;
			const E_perturbed = calculateDiscreteEnergy(
				vertices_perturbed,
				edges,
				alpha,
				beta,
				disjointPairs
			);
			differential[i][dim] = (E_perturbed - E_original) / h;
		}
	}
	console.log('Computed differential:', differential);
	return differential;
}

function calculateL2Gradient(vertices, edges, alpha, beta, disjointPairs) {
	const h = get(config).finiteDiffH; // Use config finiteDiffH
	const numVertices = vertices.length;
	const gradient = [];

	const originalEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

	for (let i = 0; i < numVertices; i++) {
		gradient[i] = [0, 0];

		for (let j = 0; j < 2; j++) {
			const originalValue = vertices[i][j];
			vertices[i][j] += h;

			const newEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

			gradient[i][j] = (newEnergy - originalEnergy) / h;

			vertices[i][j] = originalValue;
		}
	}

	return gradient;
}
