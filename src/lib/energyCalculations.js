// src/lib/energyCalculations.js
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config, subvertices as subverticesStore } from '$lib/stores';

let logging = false;

export function calculateEdgeProperties(vertices, edges) {
	const edgeLengths = [];
	const edgeTangents = [];
	const edgeMidpoints = [];

	for (const edge of edges) {
		const v1 = vertices[edge[0]];
		const v2 = vertices[edge[1]];

		// Add safety check to prevent undefined errors
		if (!v1 || !v2) {
			console.error('Invalid vertex reference:', edge, 'in vertices:', vertices);
			edgeLengths.push(0);
			edgeTangents.push([0, 0]);
			edgeMidpoints.push([0, 0]);
			continue;
		}

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

export function tangentPointKernel(p, q, T, alpha, beta) {
	// Ensure inputs are properly converted to matrices
	const p_ = math.matrix(p);
	const q_ = math.matrix(q);
	const T_ = math.matrix(T);
	const epsilon = get(config).epsilonKernel;

	// Calculate the difference vector
	const diff = math.subtract(p_, q_);
	const diffNorm = math.norm(diff) + epsilon; // Prevent division by zero
	const cross2D = T_.get([0]) * diff.get([1]) - T_.get([1]) * diff.get([0]); // 2D cross product (determinant)

	const numerator = Math.pow(Math.abs(cross2D), alpha);
	const denominator = Math.pow(diffNorm, beta);
	const result = numerator / denominator;

	// Check for NaN or Infinity
	if (!isFinite(result)) {
		console.warn('Invalid kernel result:', result, 'from inputs:', p, q, T, 'with cross2D:', cross2D, 'diffNorm:', diffNorm);
		return 0;
	}

	return result;
}

export function calculateDisjointEdgePairs(edges) {
	const numEdges = edges.length;
	const disjointPairs = [];

	for (let i = 0; i < numEdges; i++) {
		disjointPairs[i] = [];
		for (let j = 0; j < numEdges; j++) {
			if (i === j) continue;

			const edge1 = edges[i];
			const edge2 = edges[j];

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
	if (logging) {
		console.log('Calculated disjointPairs:', disjointPairs);
	}
	return disjointPairs;
}

export function calculateDiscreteKernel(vertices, edges, edgeTangents, alpha, beta, disjointPairs) {
	const numEdges = edges.length;
	const kernelMatrix = math.zeros(numEdges, numEdges);

	if (!disjointPairs || !Array.isArray(disjointPairs) || disjointPairs.length === 0) {
		console.warn('No disjoint pairs found, returning zero kernel matrix');
		return kernelMatrix;
	}

	for (let i = 0; i < numEdges; i++) {
		if (!disjointPairs[i]) {
			console.warn(`No disjoint pairs for edge ${i}`);
			continue;
		}

		for (const j of disjointPairs[i]) {
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
			if (i < edges.length && j < edges.length) {
				const kernelValue = kernelMatrix.get([i, j]);
					totalEnergy += kernelValue * edgeLengths[i] * edgeLengths[j];
			}
		}
	}
	return totalEnergy / 2; // Divide by 2 because of symmetry
}

/**
 * Prepare combined vertices and edges that include subvertices for energy calculation
 * @param {Array} mainVertices - Array of main vertex positions
 * @param {Array} mainEdges - Array of main edges
 * @param {Array} subvs - Array of subvertices
 * @returns {Object} - Combined vertices, edges, and mapping information
 */
function prepareCombinedVerticesAndEdges(mainVertices, mainEdges, subvs) {
    // Create a copy of the main vertices
    const combinedVertices = [...mainVertices.map(v => [...v])];
    const mainVertexCount = mainVertices.length;
    
    // Add subvertices to combined array first
    const subVertexIndices = new Map(); // Maps subvertex to its index in combinedVertices
    for (const sv of subvs) {
        // Add subvertex position to combined vertices
        const index = combinedVertices.length;
        combinedVertices.push([...sv.position]); // Clone the position
        
        // Store subvertex index with a key that identifies both the subvertex and its edge
        const key = `${sv.edge[0]}-${sv.edge[1]}-${sv.position[0]}-${sv.position[1]}`;
        subVertexIndices.set(key, index);
    }
    
    // Create a mapping from edge to its subvertices, sorted by position along edge
    const edgeToSubvertices = new Map();
    for (const sv of subvs) {
        const edgeKey = `${sv.edge[0]}-${sv.edge[1]}`;
        if (!edgeToSubvertices.has(edgeKey)) {
            edgeToSubvertices.set(edgeKey, []);
        }
        
        const key = `${sv.edge[0]}-${sv.edge[1]}-${sv.position[0]}-${sv.position[1]}`;
        const index = subVertexIndices.get(key);
        if (index !== undefined) {
            edgeToSubvertices.get(edgeKey).push({
                index,
                position: sv.position,
                key
            });
        }
    }
    
    // Sort subvertices along each edge
    for (const [edgeKey, svList] of edgeToSubvertices.entries()) {
        const [v1Idx, v2Idx] = edgeKey.split('-').map(Number);
        const v1 = mainVertices[v1Idx];
        
        // Sort subvertices by distance from v1
        svList.sort((a, b) => {
            const distA = Math.sqrt(Math.pow(a.position[0] - v1[0], 2) + Math.pow(a.position[1] - v1[1], 2));
            const distB = Math.sqrt(Math.pow(b.position[0] - v1[0], 2) + Math.pow(b.position[1] - v1[1], 2));
            return distA - distB;
        });
    }
    
    // Create new edge array including subvertices
    const combinedEdges = [];
    const originalEdgeIndices = new Map(); // Maps combined edge index to original edge index
    
    for (let i = 0; i < mainEdges.length; i++) {
        const edge = mainEdges[i];
        const edgeKey = `${edge[0]}-${edge[1]}`;
        const subvList = edgeToSubvertices.get(edgeKey);
        
        if (!subvList || subvList.length === 0) {
            // If no subvertices on this edge, keep the original edge
            combinedEdges.push([edge[0], edge[1]]);
            originalEdgeIndices.set(combinedEdges.length - 1, i);
        } else {
            // Create edges connecting v1 -> subv1 -> subv2 -> ... -> v2
            // First edge: v1 to first subvertex
            combinedEdges.push([edge[0], subvList[0].index]);
            originalEdgeIndices.set(combinedEdges.length - 1, i);
            
            // Middle edges: between subvertices
            for (let j = 0; j < subvList.length - 1; j++) {
                combinedEdges.push([subvList[j].index, subvList[j+1].index]);
                originalEdgeIndices.set(combinedEdges.length - 1, i);
            }
            
            // Last edge: last subvertex to v2
            combinedEdges.push([subvList[subvList.length - 1].index, edge[1]]);
            originalEdgeIndices.set(combinedEdges.length - 1, i);
        }
    }

    if (logging) {
        console.log('Combined vertices:', combinedVertices);
        console.log('Combined edges:', combinedEdges);
    }
    
    return {
        combinedVertices,
        combinedEdges,
        mainVertexCount,
        originalEdgeIndices
    };
}

/**
 * Calculate discrete energy including subvertices
 */
export function calculateDiscreteEnergyWithSubvertices(vertices, edges, alpha, beta, disjointPairs) {
    const subvs = get(subverticesStore);
    
    if (!subvs || subvs.length === 0) {
        // If no subvertices, use regular energy calculation
        return calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    }
    
    // Combine vertices and edges to include subvertices
    const { 
        combinedVertices, 
        combinedEdges
    } = prepareCombinedVerticesAndEdges(vertices, edges, subvs);
    
    // Calculate properties for combined graph
    const { edgeLengths, edgeTangents } = calculateEdgeProperties(combinedVertices, combinedEdges);
    const combinedDisjointPairs = calculateDisjointEdgePairs(combinedEdges);
    
    // Calculate kernel matrix for combined graph
    const kernelMatrix = calculateDiscreteKernel(
        combinedVertices,
        combinedEdges,
        edgeTangents,
        alpha,
        beta,
        combinedDisjointPairs
    );
    
    // Calculate energy
    let totalEnergy = 0;
    const numEdges = combinedEdges.length;
    
    for (let i = 0; i < numEdges; i++) {
        for (const j of combinedDisjointPairs[i]) {
            if (i < combinedEdges.length && j < combinedEdges.length) {
                const kernelValue = kernelMatrix.get([i, j]);
                totalEnergy += kernelValue * edgeLengths[i] * edgeLengths[j];
            }
        }
    }
    
    return totalEnergy / 2; // Divide by 2 because of symmetry
}

function calculateAnalyticalDifferential(vertices, edges, alpha, beta, disjointPairs) {
   
    return 'not implemented';}

export function calculateDifferential(vertices, edges, alpha, beta, disjointPairs) {
    const method = get(config).differentialMethod;
    const useSubvertices = get(config).useSubverticesInEnergy;
    
    if (useSubvertices) {
        return calculateDifferentialWithSubvertices(vertices, edges, alpha, beta, disjointPairs);
    }
    
    if (method === 'finiteDifference') {
        return calculateDifferentialFiniteDifference(vertices, edges, alpha, beta, disjointPairs);
    } else if (method === 'analytical') {
        return calculateAnalyticalDifferential(vertices, edges, alpha, beta, disjointPairs);
    } else {
        throw new Error('Unknown method for differential calculation');
    }
}

/**
 * Calculate the differential when including subvertices in energy calculations
 * We calculate energy with subvertices but only compute gradients for main vertices
 */
function calculateDifferentialWithSubvertices(vertices, edges, alpha, beta, disjointPairs) {
    const h = get(config).finiteDiffH;
    const numVertices = vertices.length;
    const differential = [];

    // Calculate the original energy with subvertices included
    const E_original = calculateDiscreteEnergyWithSubvertices(vertices, edges, alpha, beta, disjointPairs);

    // For each main vertex, calculate partial derivatives
    for (let i = 0; i < numVertices; i++) {
        differential[i] = [0, 0];
        for (let dim = 0; dim < 2; dim++) {
            // Create perturbed vertices array
            const vertices_perturbed = vertices.map((v) => [...v]);
            vertices_perturbed[i][dim] += h;
            
            // Calculate energy with perturbed position
            const E_perturbed = calculateDiscreteEnergyWithSubvertices(
                vertices_perturbed,
                edges,
                alpha,
                beta,
                disjointPairs
            );
            
            // Compute partial derivative using finite difference
            differential[i][dim] = (E_perturbed - E_original) / h;
        }
    }
    
    if (logging) {
        console.log('Computed differential with subvertices:', differential);
    }
    return differential;
}

function calculateDifferentialFiniteDifference(vertices, edges, alpha, beta, disjointPairs) {
    // Existing finite difference implementation
    const h = get(config).finiteDiffH;
    const numVertices = vertices.length;
    const differential = [];

    const E_original = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

    for (let i = 0; i < numVertices; i++) {
        differential[i] = [0, 0];
        for (let dim = 0; dim < 2; dim++) {
            const vertices_perturbed = vertices.map((v) => [...v]);
            vertices_perturbed[i][dim] += h;
            const E_perturbed = calculateDiscreteEnergy(
                vertices_perturbed,
                edges,
                alpha,
                beta,
                disjointPairs
            );
            differential[i][dim] = (E_perturbed - E_original) / h;
        }
    }
    if (logging) {
        console.log('Computed differential:', differential);
    }
    return differential;
}

function calculateL2Gradient(vertices, edges, alpha, beta, disjointPairs) {
	const h = get(config).finiteDiffH; // Use config finiteDiffH
	const numVertices = vertices.length;
	const gradient = [];

	const originalEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

	for (let i = 0; i < numVertices; i++) {
		gradient[i] = [0, 0];

		for (let j = 0; j < 2; j++) {
			const originalValue = vertices[i][j];
			vertices[i][j] += h;

			const newEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

			gradient[i][j] = (newEnergy - originalEnergy) / h;

			vertices[i][j] = originalValue;
		}
	}

	return gradient;
}