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
	const p_ = math.matrix(p);
	const q_ = math.matrix(q);
	const T_ = math.matrix(T);

	const diff = math.subtract(p_, q_);
	const cross = math.cross(T_, diff);

	const numerator = Math.pow(math.norm(cross), alpha);
	const denominator = Math.pow(math.norm(diff), beta);
	return math.divide(numerator, denominator); // optimized using math divide
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
			const tx_cross_p_minus_q = math.cross([...ti, 0], [...p_minus_q, 0]); // 2D cross product (z-component)

			const numerator = Math.pow(Math.abs(tx_cross_p_minus_q[2]), beta); // Use abs for 2D
			const denominator = Math.pow(math.norm(p_minus_q), beta);

			const k_beta = numerator / denominator;
			const edgeEnergy = k_beta * edgeLengths[i] * edgeLengths[j];
			totalEnergy += edgeEnergy;

			// --- Gradient Calculation (for animation) ---
			const d_k_beta_d_xi = math.multiply(
				(beta * k_beta) / math.norm(p_minus_q),
				math.subtract(
					math.multiply(math.dot(ti, p_minus_q) / math.norm(p_minus_q), ti),
					math.multiply(
						(tx_cross_p_minus_q[2] * tx_cross_p_minus_q[2]) / (numerator + 1e-9),
						p_minus_q
					) //add small number to prevent NaN
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
