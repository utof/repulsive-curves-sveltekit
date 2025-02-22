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
function calculateE_adj(vertices, edges) {
	const numVertices = vertices.length;
	const E_adj = [];
	for (let i = 0; i < numVertices; i++) {
		E_adj.push([]); // Initialize an empty array for each vertex
	}

	for (let i = 0; i < edges.length; i++) {
		const edge = edges[i];
		E_adj[edge[0]].push(i); // Add edge index to both vertices' adjacency lists
		E_adj[edge[1]].push(i);
	}
	return E_adj;
}

export function computeGradient(vertices, edges, alpha, beta, disjointPairs, edgeProps) {
	const numVertices = vertices.length;
	const numEdges = edges.length;
	const gradient = [];
	for (let i = 0; i < numVertices; i++) {
		gradient.push([0, 0]); // Initialize gradient for each vertex
	}

	const E_adj = calculateE_adj(vertices, edges); // Calculate E_adj

	for (let p = 0; p < numVertices; p++) {
		let derivP = [0, 0]; // Gradient for the current vertex

		for (const I of E_adj[p]) {
			// Edges adjacent to p
			for (const J of disjointPairs[I]) {
				// Edges disjoint from I
				for (let i = 0; i < 2; i++) {
					for (let j = 0; j < 2; j++) {
						if (edges[I][i] === p) {
							const i1 = edges[I][i];
							const i2 = edges[I][(i + 1) % 2];
							const j1 = edges[J][j];

							// --- Case 1,1 ---
							const crossTerm = [
								vertices[i2][0] - vertices[j1][0], // x component
								vertices[i2][1] - vertices[j1][1] // y component
							];

							const denomDiff = [
								vertices[i1][0] - vertices[j1][0],
								vertices[i1][1] - vertices[j1][1]
							];

							const crossTerm_cross_e1 = -crossTerm[1];
							const crossTerm_cross_e2 = crossTerm[0];

							const TMat = [
								[crossTerm_cross_e1, 0],
								[0, crossTerm_cross_e2]
							];

							const term1 = math.multiply(
								(1 - alpha) * Math.pow(edgeProps.edgeLengths[I], -alpha - 1),
								[vertices[i1][0] - vertices[i2][0], vertices[i1][1] - vertices[i2][1]],
								Math.pow(math.norm([crossTerm_cross_e1, crossTerm_cross_e2]), alpha),
								Math.pow(math.norm(denomDiff), -beta)
							);

							const term2 = math.multiply(
								alpha * Math.pow(edgeProps.edgeLengths[I], 1 - alpha),
								Math.pow(math.norm(denomDiff), -beta),
								Math.pow(math.norm([crossTerm_cross_e1, crossTerm_cross_e2]), alpha - 2),
								math.multiply([crossTerm_cross_e1, crossTerm_cross_e2], TMat)
							);

							const term3 = math.multiply(
								-beta * Math.pow(edgeProps.edgeLengths[I], 1 - alpha),
								Math.pow(math.norm(denomDiff), -beta - 2),
								Math.pow(math.norm([crossTerm_cross_e1, crossTerm_cross_e2]), alpha),
								denomDiff
							);

							const termSum11 = math.multiply(
								0.25 * edgeProps.edgeLengths[J],
								math.add(math.add(term1, term2), term3)
							);
							derivP = math.add(derivP, termSum11);

							// --- Case 2,1 ---
							const denomDiff21 = [
								vertices[i2][0] - vertices[j1][0],
								vertices[i2][1] - vertices[j1][1]
							];

							const term1_21 = math.multiply(
								(1 - alpha) * Math.pow(edgeProps.edgeLengths[I], -alpha - 1),
								[vertices[i1][0] - vertices[i2][0], vertices[i1][1] - vertices[i2][1]],
								Math.pow(math.norm([crossTerm_cross_e1, crossTerm_cross_e2]), alpha),
								Math.pow(math.norm(denomDiff21), -beta)
							);

							const term2_21 = math.multiply(
								alpha * Math.pow(edgeProps.edgeLengths[I], 1 - alpha),
								Math.pow(math.norm(denomDiff21), -beta),
								Math.pow(math.norm([crossTerm_cross_e1, crossTerm_cross_e2]), alpha - 2),
								math.multiply([crossTerm_cross_e1, crossTerm_cross_e2], TMat)
							);

							const termSum21 = math.multiply(
								0.25 * edgeProps.edgeLengths[J],
								math.add(term1_21, term2_21)
							);
							derivP = math.add(derivP, termSum21);

							// --- Case J: 1,1 ---
							const TJ = edgeProps.edgeTangents[J];
							const TMatJ = [
								[-TJ[1], 0],
								[0, TJ[0]]
							];

							const denomDiffJ11 = [
								vertices[i1][0] - vertices[j1][0],
								vertices[i1][1] - vertices[j1][1]
							];

							const crossTermJ11 = [
								-TJ[1] * denomDiffJ11[1], // -TJy * diffy
								TJ[0] * denomDiffJ11[0] //  TJX * diffx
							];

							const term1J11 = math.multiply(
								Math.pow(math.norm(crossTermJ11), alpha),
								Math.pow(math.norm(denomDiffJ11), -beta),
								denomDiffJ11,
								1 / edgeProps.edgeLengths[I]
							);

							const term2J11 = math.multiply(
								math.multiply(crossTermJ11, TMatJ),
								alpha * edgeProps.edgeLengths[I],
								Math.pow(math.norm(denomDiffJ11), -beta),
								Math.pow(math.norm(crossTermJ11), alpha - 2)
							);

							const term3J11 = math.multiply(
								-beta * edgeProps.edgeLengths[I],
								Math.pow(math.norm(denomDiffJ11), -beta - 2),
								Math.pow(math.norm(crossTermJ11), alpha),
								denomDiffJ11
							);

							const termSumJ11 = math.multiply(
								0.25 * edgeProps.edgeLengths[J],
								math.add(math.add(term1J11, term2J11), term3J11)
							);
							derivP = math.add(derivP, termSumJ11);

							// --- Case J: 1,2 ---
							const denomDiffJ12 = [
								vertices[i2][0] - vertices[j1][0],
								vertices[i2][1] - vertices[j1][1]
							];

							const crossTermJ12 = [
								-TJ[1] * denomDiffJ12[1], // -TJy * diffy
								TJ[0] * denomDiffJ12[0] //  TJX * diffx
							];

							const termJ12 = math.multiply(
								0.25 * (edgeProps.edgeLengths[J] / edgeProps.edgeLengths[I]),
								Math.pow(math.norm(crossTermJ12), alpha),
								Math.pow(math.norm(denomDiffJ12), -beta),
								[vertices[i1][0] - vertices[i2][0], vertices[i1][1] - vertices[i2][1]]
							);

							derivP = math.add(derivP, termJ12);
						}
					}
				}
			}
		}
		gradient[p][0] = derivP[0];
		gradient[p][1] = derivP[1];
	}
	return gradient;
}
