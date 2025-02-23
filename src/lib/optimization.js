// src/lib/optimization.js
import { calculateL2Gradient, calculateEdgeProperties } from '$lib/energyCalculations';
import { computePreconditionedGradient } from '$lib/innerProduct';

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
	console.log('Gradient Descent Step - Initial vertices:', vertices);
	const { edgeTangents } = calculateEdgeProperties(vertices, edges);
	const L2Gradient = calculateL2Gradient(vertices, edges, alpha, beta, disjointPairs);
	const sigma = beta / alpha - 1.5;
	console.log('Sigma:', sigma);
	const gradient = computePreconditionedGradient(vertices, edges, edgeTangents, sigma, L2Gradient);

	const newVertices = vertices.map((vertex, i) => {
		const gradX = gradient[i][0];
		const gradY = gradient[i][1];
		const newX = Math.max(0, Math.min(width, vertex[0] - stepSize * (isNaN(gradX) ? 0 : gradX)));
		const newY = Math.max(0, Math.min(height, vertex[1] - stepSize * (isNaN(gradY) ? 0 : gradY)));
		console.log(
			`Vertex[${i}] update: [${vertex[0]}, ${vertex[1]}] -> [${newX}, ${newY}] with grad [${gradX}, ${gradY}]`
		);
		return [newX, newY];
	});

	return newVertices;
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
			console.log(`Optimization Step ${currentIteration}`);
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
