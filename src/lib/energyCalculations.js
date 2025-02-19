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
	// 2D cross product is just T.x * diff.y - T.y * diff.x
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
			if (
				i === j ||
				edges[i][0] === edges[j][0] ||
				edges[i][0] === edges[j][1] ||
				edges[i][1] === edges[j][0] ||
				edges[i][1] === edges[j][1]
			) {
				continue; // Skip neighboring edges
			}

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
		}
	}
	return kernelMatrix;
}

export function discreteTangentPointEnergy(vertices, edges, alpha, beta) {
	let totalEnergy = 0;
	const vertexData = vertices.map(() => ({ energy: 0, count: 0, gradient: [0, 0] }));

	const edgeLengths = edges.map((edge) => math.distance(vertices[edge[0]], vertices[edge[1]]));
	const edgeTangents = edges.map((edge, i) =>
		math.divide(math.subtract(vertices[edge[1]], vertices[edge[0]]), edgeLengths[i])
	);

	for (let i = 0; i < edges.length; i++) {
		for (let j = 0; j < edges.length; j++) {
			if (i === j || edges[i].some((v) => edges[j].includes(v))) continue;

			const ti = edgeTangents[i];
			const xi = vertices[edges[i][0]];
			const xj = vertices[edges[j][0]];

			const p_minus_q = math.subtract(xi, xj);
			// 2D cross product is just T.x * diff.y - T.y * diff.x
			const cross2D = ti[0] * p_minus_q[1] - ti[1] * p_minus_q[0];

			const numerator = Math.pow(Math.abs(cross2D), beta);
			const denominator = Math.pow(math.norm(p_minus_q), beta);

			const k_beta = numerator / denominator;
			const edgeEnergy = k_beta * edgeLengths[i] * edgeLengths[j];
			totalEnergy += edgeEnergy;

			// --- Gradient Calculation (for animation) ---
			const d_k_beta_d_xi = math.multiply(
				(beta * k_beta) / math.norm(p_minus_q),
				math.subtract(
					math.multiply(math.dot(ti, p_minus_q) / math.norm(p_minus_q), ti),
					math.multiply((cross2D * cross2D) / (numerator + 1e-9), p_minus_q) //add small number to prevent NaN
				)
			);

			// Accumulate gradient and energy per vertex
			vertexData[edges[i][0]].gradient = math.add(
				vertexData[edges[i][0]].gradient,
				math.multiply(edgeLengths[i] * edgeLengths[j], d_k_beta_d_xi)
			);
			vertexData[edges[i][1]].gradient = math.subtract(
				vertexData[edges[i][1]].gradient,
				math.multiply(edgeLengths[i] * edgeLengths[j], d_k_beta_d_xi)
			);

			vertexData[edges[i][0]].energy += edgeEnergy;
			vertexData[edges[i][0]].count += 1;
			vertexData[edges[i][1]].energy += edgeEnergy;
			vertexData[edges[i][1]].count += 1;
		}
	}

	for (let i = 0; i < vertexData.length; i++) {
		if (vertexData[i].count > 0) {
			vertexData[i].energy /= vertexData[i].count;
		}
	}

	return { totalEnergy, vertexData };
}
