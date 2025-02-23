// src/lib/optimization.js
import { calculateL2Gradient } from '$lib/energyCalculations';

export function gradientDescentStep(
	vertices,
	edges,
	alpha,
	beta,
	disjointPairs,
	stepSize,
	width,
	height
) {
	const gradient = calculateL2Gradient(vertices, edges, alpha, beta, disjointPairs);

	return vertices.map((vertex, i) => {
		const newX = Math.max(0, Math.min(width, vertex[0] - stepSize * gradient[i][0]));
		const newY = Math.max(0, Math.min(height, vertex[1] - stepSize * gradient[i][1]));
		return [newX, newY];
	});
}

export function createOptimizer(
	vertices,
	edges,
	alpha,
	beta,
	disjointPairs,
	stepSize,
	width,
	height,
	maxIterations,
	onUpdate
) {
	let currentIteration = 0;
	let intervalId = null;

	const step = () => {
		if (currentIteration < maxIterations) {
			const newVertices = gradientDescentStep(
				vertices,
				edges,
				alpha,
				beta,
				disjointPairs,
				stepSize,
				width,
				height
			);
			for (let i = 0; i < vertices.length; i++) {
				vertices[i][0] = newVertices[i][0];
				vertices[i][1] = newVertices[i][1];
			}
			currentIteration++;
			onUpdate();
		} else {
			stop();
		}
	};

	const start = () => {
		currentIteration = 0;
		if (!intervalId) {
			intervalId = setInterval(step, 20);
		}
	};

	const stop = () => {
		if (intervalId) {
			clearInterval(intervalId);
			intervalId = null;
		}
	};

	return { start, stop, step };
}
