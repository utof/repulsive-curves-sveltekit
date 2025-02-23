// src/lib/energyCalculations.js
import * as math from 'mathjs';

let logging = false;

export function calculateEdgeProperties(vertices, edges) {
	const edgeLengths = [];
	const edgeTangents = [];
	const edgeMidpoints = [];

	for (const edge of edges) {
		const v1 = vertices[edge[0]];
		const v2 = vertices[edge[1]];

		const dx = v2[0] - v1[0];
		const dy = v2[1] - v1[1];
		const length = Math.sqrt(dx * dx + dy * dy);
		edgeLengths.push(length);

		const unitTangent = length > 0 ? [dx / length, dy / length] : [0, 0];
		edgeTangents.push(unitTangent);

		const midpoint = [
			isNaN(v1[0]) || isNaN(v2[0]) ? 0 : (v1[0] + v2[0]) / 2,
			isNaN(v1[1]) || isNaN(v2[1]) ? 0 : (v1[1] + v2[1]) / 2
		];
		edgeMidpoints.push(midpoint);

		if (logging) {
			console.log(
				`Edge [${edge[0]}, ${edge[1]}]: length = ${length}, tangent = ${unitTangent}, midpoint = ${midpoint}`
			);
		}
	}

	if (logging) {
		console.log('Edge lengths:', edgeLengths);
		console.log('Unit tangents:', edgeTangents);
		console.log('Midpoints:', edgeMidpoints);
	}

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

	// console.log('Kernel calc:', {
	// 	p: p_.toArray(),
	// 	q: q_.toArray(),
	// 	T: T_.toArray(),
	// 	cross2D,
	// 	diffNorm,
	// 	numerator,
	// 	denominator,
	// 	result
	// });

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

// Function to calculate the differential (not the L2 gradient)
export function calculateDifferential(vertices, edges, alpha, beta, disjointPairs) {
	const numVertices = vertices.length;
	const differential = [];

	for (let i = 0; i < numVertices; i++) {
		differential[i] = [0, 0]; // Initialize for x and y components
	}

	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);

	for (let I = 0; I < edges.length; I++) {
		for (const J of disjointPairs[I]) {
			const combinations = [
				[vertices[edges[I][0]], vertices[edges[J][0]]],
				[vertices[edges[I][0]], vertices[edges[J][1]]],
				[vertices[edges[I][1]], vertices[edges[J][0]]],
				[vertices[edges[I][1]], vertices[edges[J][1]]]
			];

			for (let k = 0; k < combinations.length; k++) {
				const [p, q] = combinations[k];
				const T = edgeTangents[I];

				// Compute partial derivatives of the kernel
				const h = 1e-6;
				const kernelValue = tangentPointKernel(p, q, T, alpha, beta);

				// Partial derivatives with respect to p_x, p_y, q_x, q_y
				const dp_x = (tangentPointKernel([p[0] + h, p[1]], q, T, alpha, beta) - kernelValue) / h;
				const dp_y = (tangentPointKernel([p[0], p[1] + h], q, T, alpha, beta) - kernelValue) / h;
				const dq_x = (tangentPointKernel(p, [q[0] + h, q[1]], T, alpha, beta) - kernelValue) / h;
				const dq_y = (tangentPointKernel(p, [q[0], q[1] + h], T, alpha, beta) - kernelValue) / h;

				// Determine which vertices these partials correspond to
				let p_vertexIndex, q_vertexIndex;
				if (k < 2) {
					// p is from edge I's first vertex
					p_vertexIndex = edges[I][0];
				} else {
					// p is from edge I's second vertex
					p_vertexIndex = edges[I][1];
				}
				if (k % 2 === 0) {
					// q is from edge J's first vertex
					q_vertexIndex = edges[J][0];
				} else {
					// q is from edge J's second vertex
					q_vertexIndex = edges[J][1];
				}

				// Accumulate the contributions (scaled by edge lengths)
				differential[p_vertexIndex][0] += (dp_x * edgeLengths[I] * edgeLengths[J]) / 4;
				differential[p_vertexIndex][1] += (dp_y * edgeLengths[I] * edgeLengths[J]) / 4;
				differential[q_vertexIndex][0] += (dq_x * edgeLengths[I] * edgeLengths[J]) / 4;
				differential[q_vertexIndex][1] += (dq_y * edgeLengths[I] * edgeLengths[J]) / 4;
			}
		}
	}

	return differential;
}
// for later use.
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
