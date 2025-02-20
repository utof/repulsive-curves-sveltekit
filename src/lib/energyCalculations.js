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
	const diffNorm = math.norm(diff);
	// if (diffNorm < 1e-10) return 0; // Points too close together

	// cross = T x (p - q)
	const cross2D = T_.get([0]) * diff.get([1]) - T_.get([1]) * diff.get([0]);

	// Scale the result for better visibility
	const epsilon = 1e-6; // Prevent division by zero
	const numerator = Math.pow(Math.abs(cross2D), alpha);
	const denominator = Math.pow(diffNorm + epsilon, beta);
	const result = numerator / denominator; // Scale by 100 for better visibility

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

	return result;
}

export function calculateDisjointEdgePairs(edges) {
	const numEdges = edges.length;
	const disjointPairs = [];

	for (let i = 0; i < numEdges; i++) {
		disjointPairs[i] = []; // Initialize the array for edge i
		for (let j = 0; j < numEdges; j++) {
			if (i === j) continue; // Don't compare an edge to itself

			const edge1 = edges[i];
			const edge2 = edges[j];

			// Check if the edges share any vertices.  If they don't share any,
			// they are disjoint.
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
	console.log('Calculated disjointPairs:', disjointPairs); // Add this
	return disjointPairs;
}

export function calculateDiscreteKernel(vertices, edges, edgeTangents, alpha, beta, disjointPairs) {
	const numEdges = edges.length;
	const kernelMatrix = math.zeros(numEdges, numEdges);

	// Add defensive checks for disjointPairs
	if (!disjointPairs || !Array.isArray(disjointPairs) || disjointPairs.length === 0) {
		console.warn('No disjoint pairs found, returning zero kernel matrix');
		return kernelMatrix;
	}

	for (let i = 0; i < numEdges; i++) {
		// Add defensive check for disjointPairs[i]
		if (!disjointPairs[i]) {
			console.warn(`No disjoint pairs for edge ${i}`);
			continue;
		}

		for (const j of disjointPairs[i]) {
			// Defensive check: Ensure i and j are valid indices
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
			// No need for the neighbor check anymore!
			if (i < edges.length && j < edges.length) {
				// Add defensive check
				const kernelValue = kernelMatrix.get([i, j]);
				totalEnergy += kernelValue * edgeLengths[i] * edgeLengths[j];
			}
		}
	}
	return totalEnergy / 2; // Divide by 2 because of symmetry
}
