// src/lib/optimization.js
import {
	calculateDifferential,
	calculateEdgeProperties,
	calculateDiscreteEnergy
} from '$lib/energyCalculations'; // Corrected import
import { computePreconditionedGradient } from '$lib/innerProduct';

import * as math from 'mathjs';

export function gradientDescentStep(
	vertices,
	edges,
	alpha,
	beta,
	disjointPairs,
	a_const = 0.1,
	b_const = 0.5,
	max_line_search = 20
) {
	console.log('Gradient Descent Step - Initial vertices:', JSON.stringify(vertices));
	const { edgeTangents } = calculateEdgeProperties(vertices, edges);

	const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
	console.log('Differential:', differential);

	const gradient = computePreconditionedGradient(
		vertices,
		edges,
		edgeTangents,
		alpha,
		beta,
		differential
	);
	console.log('Preconditioned Gradient:', gradient);

	const d = gradient.map((g, i) => (i === 0 ? [0, 0] : [-g[0], -g[1]]));
	const d_norm = Math.sqrt(d.flat().reduce((sum, val) => sum + val * val, 0)) || 1;
	const d_normalized = d.map(([dx, dy]) => [dx / d_norm, dy / d_norm]);
	console.log('Normalized descent direction:', d_normalized);

	const differentialFlat = differential.flat();
	const dFlat = d_normalized.flat();
	const slope = math.dot(differentialFlat, dFlat);
	console.log('Slope (nabla E^T d):', slope);

	const E_old = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
	console.log('Initial energy:', E_old);

	let t = 1.0;
	for (let i = 0; i < max_line_search; i++) {
		const vertices_new = vertices.map((vertex, i) => {
			if (i === 0) return vertex;
			return [vertex[0] + t * d_normalized[i][0], vertex[1] + t * d_normalized[i][1]];
		});

		const E_new = calculateDiscreteEnergy(vertices_new, edges, alpha, beta, disjointPairs);
		console.log(`Line search iteration ${i}, t=${t}, E_new=${E_new}`);

		if (E_new <= E_old + a_const * t * slope) {
			console.log('Line search converged at t=', t);
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
	onUpdate
) {
	if (typeof onUpdate !== 'function') {
		throw new Error('onUpdate must be a function');
	}

	let currentIteration = 0;
	let intervalId = null;

	const optimizer = {
		step: () => {
			if (currentIteration < maxIterations) {
				console.log(`Optimization Step ${currentIteration} - Starting`);
				const newVertices = gradientDescentStep(vertices, edges, alpha, beta, disjointPairs);
				vertices.forEach((v, i) => {
					v[0] = newVertices[i][0];
					v[1] = newVertices[i][1];
				});
				currentIteration++;
				console.log(`Step ${currentIteration} - Vertices updated:`, JSON.stringify(vertices));
				onUpdate();
				console.log(`Step ${currentIteration} - onUpdate called`);
			} else {
				optimizer.stop();
				console.log('Optimization stopped - Max iterations reached');
			}
		},

		start: () => {
			currentIteration = 0;
			if (!intervalId) {
				console.log('Starting optimization');
				intervalId = setInterval(optimizer.step, 20);
			}
		},

		stop: () => {
			if (intervalId) {
				clearInterval(intervalId);
				intervalId = null;
				console.log('Optimization stopped');
			}
		}
	};

	return optimizer;
}
