// src/lib/optimization.js
import {
	calculateDifferential,
	calculateEdgeProperties,
	calculateDiscreteEnergy
} from '$lib/energyCalculations';
import { computePreconditionedGradient } from '$lib/innerProduct';
import * as math from 'mathjs';

function projectConstraints(
	vertices,
	edges,
	initialEdgeLengths,
	maxIterations = 10,
	tolerance = 1e-7
) {
	for (let iter = 0; iter < maxIterations; iter++) {
		let maxError = 0;
		for (let I = 0; I < edges.length; I++) {
			const [i, j] = edges[I];
			const vi = vertices[i];
			const vj = vertices[j];
			const d = [vj[0] - vi[0], vj[1] - vi[1]];
			const currentLength = Math.sqrt(d[0] * d[0] + d[1] * d[1]);
			const targetLength = initialEdgeLengths[I];
			const error = currentLength - targetLength;
			if (Math.abs(error) > tolerance) {
				const correction = (error / currentLength) * 0.5;
				const delta = [d[0] * correction, d[1] * correction];
				vertices[i][0] += delta[0];
				vertices[i][1] += delta[1];
				vertices[j][0] -= delta[0];
				vertices[j][1] -= delta[1];
				maxError = Math.max(maxError, Math.abs(error));
			}
		}
		if (maxError < tolerance) break;
	}
}

export function gradientDescentStep(
	vertices,
	edges,
	alpha,
	beta,
	disjointPairs,
	initialEdgeLengths,
	a_const = 0.1,
	b_const = 0.5,
	max_line_search = 20
) {
	console.log('Gradient Descent Step - Initial vertices:', JSON.stringify(vertices));
	const { edgeTangents } = calculateEdgeProperties(vertices, edges);

	const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
	const gradient = computePreconditionedGradient(
		vertices,
		edges,
		edgeTangents,
		alpha,
		beta,
		differential
	);

	// Allow all vertices to move (remove fixed vertex)
	const d = gradient.map(([gx, gy]) => [-gx, -gy]);
	const d_norm = Math.sqrt(d.flat().reduce((sum, val) => sum + val * val, 0)) || 1;
	const d_normalized = d.map(([dx, dy]) => [dx / d_norm, dy / d_norm]);

	const differentialFlat = differential.flat();
	const dFlat = d_normalized.flat();
	const slope = math.dot(differentialFlat, dFlat);

	const E_old = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

	let t = 1.0;
	for (let i = 0; i < max_line_search; i++) {
		const vertices_new = vertices.map((vertex, idx) => [
			vertex[0] + t * d_normalized[idx][0],
			vertex[1] + t * d_normalized[idx][1]
		]);

		// Project onto edge length constraints
		// projectConstraints(vertices_new, edges, initialEdgeLengths);

		const E_new = calculateDiscreteEnergy(vertices_new, edges, alpha, beta, disjointPairs);
		if (E_new <= E_old + a_const * t * slope) {
			return vertices_new;
		}
		t *= b_const;
	}

	console.warn('Line search did not converge');
	return vertices;
}

export function createOptimizer(
	vertices,
	edges,
	alpha,
	beta,
	disjointPairs,
	maxIterations,
	onUpdate,
	initialEdgeLengths
) {
	if (typeof onUpdate !== 'function') throw new Error('onUpdate must be a function');

	let currentIteration = 0;
	let intervalId = null;

	const optimizer = {
		step: () => {
			if (currentIteration < maxIterations) {
				const newVertices = gradientDescentStep(
					vertices,
					edges,
					alpha,
					beta,
					disjointPairs,
					initialEdgeLengths
				);
				vertices.forEach((v, i) => {
					v[0] = newVertices[i][0];
					v[1] = newVertices[i][1];
				});
				currentIteration++;
				onUpdate();
			} else {
				optimizer.stop();
			}
		},
		start: () => {
			currentIteration = 0;
			if (!intervalId) intervalId = setInterval(optimizer.step, 20);
		},
		stop: () => {
			if (intervalId) {
				clearInterval(intervalId);
				intervalId = null;
			}
		}
	};

	return optimizer;
}
