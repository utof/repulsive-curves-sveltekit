import * as math from 'mathjs';

export function calculateEdgeProperties(vertices, edges) {
	const edgeLengths = [];
	const edgeTangents = [];
	const edgeMidpoints = [];

	for (const edge of edges) {
		const v1 = vertices[edge[0]];
		const v2 = vertices[edge[1]];

		// Calculate edge length (l)
		const dx = v2[0] - v1[0];
		const dy = v2[1] - v1[1];
		const length = Math.sqrt(dx * dx + dy * dy);
		edgeLengths.push(length);

		// Calculate unit tangent (T)
		const unitTangent = [dx / length, dy / length];
		edgeTangents.push(unitTangent);

		// Calculate midpoint
		const midpoint = [(v1[0] + v2[0]) / 2, (v1[1] + v2[1]) / 2];
		edgeMidpoints.push(midpoint);
	}

	console.log('Edge lengths:', edgeLengths);
	console.log('Unit tangents:', edgeTangents);
	console.log('Midpoints:', edgeMidpoints);

	return { edgeLengths, edgeTangents, edgeMidpoints };
}

function tangentPointKernel(p, q, T, alpha, beta) {
	// For 2D vectors, cross product is just determinant of 2x2 matrix
	const p_ = math.matrix(p);
	const q_ = math.matrix(q);
	const T_ = math.matrix(T);

	const diff = math.subtract(p_, q_);
	// cross = T x (p - q)
	const cross2D = T_.get([0]) * diff.get([1]) - T_.get([1]) * diff.get([0]);

	const numerator = Math.pow(Math.abs(cross2D), alpha);
	const denominator = Math.pow(math.norm(diff), beta);
	return numerator / denominator;
}

export function calculateDiscreteKernel(vertices, edges, edgeTangents, alpha, beta) {
	const numEdges = edges.length;
	const kernelMatrix = math.zeros(numEdges, numEdges);

	for (let i = 0; i < numEdges; i++) {
		for (let j = 0; j < numEdges; j++) {
			let sum = 0;
			const combinations = [
				[vertices[edges[i][0]], vertices[edges[j][0]]],
				[vertices[edges[i][0]], vertices[edges[j][1]]],
				[vertices[edges[i][1]], vertices[edges[j][0]]],
				[vertices[edges[i][1]], vertices[edges[j][1]]]
			];
			// TODO optimize because matrix is symmetric.
			for (const [p, q] of combinations) {
				sum += tangentPointKernel(p, q, edgeTangents[i], alpha, beta);
			}

			kernelMatrix.set([i, j], sum / 4);
		}
	}
	return kernelMatrix;
}

export function calculateDiscreteEnergy(vertices, edges, alpha, beta) {
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
	const kernelMatrix = calculateDiscreteKernel(vertices, edges, edgeTangents, alpha, beta);

	let totalEnergy = 0;
	const numEdges = edges.length;

	for (let i = 0; i < numEdges; i++) {
		for (let j = 0; j < numEdges; j++) {
			// Ensure we're only considering non-neighboring edges, as per the paper.
			if (
				i === j ||
				edges[i][0] === edges[j][0] ||
				edges[i][0] === edges[j][1] ||
				edges[i][1] === edges[j][0] ||
				edges[i][1] === edges[j][1]
			) {
				continue;
			}

			const kernelValue = kernelMatrix.get([i, j]);
			totalEnergy += kernelValue * edgeLengths[i] * edgeLengths[j];
		}
	}
	return totalEnergy;
}
