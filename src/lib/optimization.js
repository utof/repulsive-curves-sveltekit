import { calculateDiscreteEnergy } from '$lib/energyCalculations'; // Import from energyCalculations

function calculateL2Gradient(vertices, edges, alpha, beta, disjointPairs) {
	const h = 0.0001; // Small change for finite differences
	const numVertices = vertices.length;
	const gradient = [];

	const originalEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

	for (let i = 0; i < numVertices; i++) {
		gradient[i] = [0, 0]; // Initialize gradient for this vertex

		for (let j = 0; j < 2; j++) {
			// x and y coordinates
			// Perturb the vertex coordinate
			const originalValue = vertices[i][j];
			vertices[i][j] += h;

			// Recalculate the energy
			const newEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

			// Approximate the partial derivative
			gradient[i][j] = (newEnergy - originalEnergy) / h;

			// Restore the original value
			vertices[i][j] = originalValue;
		}
	}

	return gradient;
}

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

	const newVertices = vertices.map((vertex, i) => {
		const newX = Math.max(0, Math.min(width, vertex[0] - stepSize * gradient[i][0]));
		const newY = Math.max(0, Math.min(height, vertex[1] - stepSize * gradient[i][1]));
		return [newX, newY];
	});

	return newVertices; // Return the *new* vertex positions
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
			// Update vertices in place (important for reactivity)
			for (let i = 0; i < vertices.length; i++) {
				vertices[i][0] = newVertices[i][0];
				vertices[i][1] = newVertices[i][1];
			}

			currentIteration++;
			onUpdate(); // Call the callback
		} else {
			stop();
		}
	};

	const start = () => {
		currentIteration = 0;
		if (!intervalId) {
			intervalId = setInterval(step, 20); // Adjust interval for speed
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
