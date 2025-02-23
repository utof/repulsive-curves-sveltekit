// src/lib/optimization.js
import { calculateL2Gradient, calculateEdgeProperties } from '$lib/energyCalculations';
import { computePreconditionedGradient } from '$lib/innerProduct';

export function gradientDescentStep(vertices, edges, alpha, beta, disjointPairs, stepSize) {
	console.log('Gradient Descent Step - Initial vertices:', JSON.stringify(vertices));
	const { edgeTangents } = calculateEdgeProperties(vertices, edges);
	const L2Gradient = calculateL2Gradient(vertices, edges, alpha, beta, disjointPairs);
	console.log('L2 Gradient:', JSON.stringify(L2Gradient));
	const gradient = computePreconditionedGradient(
		vertices,
		edges,
		edgeTangents,
		alpha,
		beta,
		L2Gradient
	);
	console.log('Preconditioned Gradient:', JSON.stringify(gradient));

	const newVertices = vertices.map((vertex, i) => {
		if (i === 0) {
			console.log(`Vertex[${i}] fixed: [${vertex[0]}, ${vertex[1]}]`);
			return [vertex[0], vertex[1]];
		}
		const gradX = gradient[i][0];
		const gradY = gradient[i][1];
		const newX = vertex[0] - stepSize * (isNaN(gradX) ? 0 : gradX);
		const newY = vertex[1] - stepSize * (isNaN(gradY) ? 0 : gradY);
		console.log(
			`Vertex[${i}] update: [${vertex[0]}, ${vertex[1]}] -> [${newX}, ${newY}] with grad [${gradX}, ${gradY}]`
		);
		return [newX, newY];
	});

	console.log('New vertices:', JSON.stringify(newVertices));
	return newVertices;
}

export function createOptimizer(
	vertices,
	edges,
	alpha,
	beta,
	disjointPairs,
	stepSize,
	maxIterations,
	onUpdate
) {
	let currentIteration = 0;
	let intervalId = null;

	const step = () => {
		if (currentIteration < maxIterations) {
			console.log(`Optimization Step ${currentIteration} - Starting`);
			const newVertices = gradientDescentStep(
				vertices,
				edges,
				alpha,
				beta,
				disjointPairs,
				stepSize
			);
			console.log(`Step ${currentIteration} - New vertices computed:`, JSON.stringify(newVertices));
			for (let i = 0; i < vertices.length; i++) {
				vertices[i][0] = newVertices[i][0];
				vertices[i][1] = newVertices[i][1];
			}
			currentIteration++;
			console.log(`Step ${currentIteration} - Vertices updated:`, JSON.stringify(vertices));
			onUpdate();
			console.log(`Step ${currentIteration} - onUpdate called`);
		} else {
			stop();
			console.log('Optimization stopped - Max iterations reached');
		}
	};

	const start = () => {
		currentIteration = 0;
		if (!intervalId) {
			console.log('Starting optimization');
			intervalId = setInterval(step, 20);
		}
	};

	const stop = () => {
		if (intervalId) {
			clearInterval(intervalId);
			intervalId = null;
			console.log('Optimization stopped');
		}
	};

	return { start, stop, step };
}
