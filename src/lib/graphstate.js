// src/lib/graphState.js
import {
	calculateEdgeProperties,
	calculateDisjointEdgePairs,
	calculateDiscreteKernel,
	calculateDiscreteEnergy,
	calculateDiscreteEnergyWithSubvertices
} from '$lib/energyCalculations';
import { generateSubvertices } from '$lib/graphUtils';
import { subvertices, config } from '$lib/stores';
import { get } from 'svelte/store';

/**
 * Initialize the kernel state with vertices and edges
 * Supports both 2D and 3D based on the config dimension
 */
export function initializeKernelState(vertices, edges, alpha, beta) {
	const disjointPairs = calculateDisjointEdgePairs(edges);
	const edgeProps = calculateEdgeProperties(vertices, edges);
	const kernelMatrix = calculateDiscreteKernel(
		vertices,
		edges,
		edgeProps.edgeTangents,
		alpha,
		beta,
		disjointPairs
	);
	
	// Calculate energy based on configuration
	const useSubvertices = get(config).useSubverticesInEnergy;
	const discreteEnergy = useSubvertices
		? calculateDiscreteEnergyWithSubvertices(vertices, edges, alpha, beta, disjointPairs)
		: calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
	
	// Generate subvertices based on current dimension
	const newSubvertices = generateSubvertices(vertices, edges);
	subvertices.set(newSubvertices);

	return {
		kernelMatrix,
		discreteEnergy,
		edgeProps,
		disjointPairs
	};
}

/**
 * Update the kernel state after vertex positions have changed
 * Works correctly with both 2D and 3D
 */
export function updateKernelState(vertices, edges, alpha, beta, disjointPairs) {
	const edgeProps = calculateEdgeProperties(vertices, edges);
	const kernelMatrix = calculateDiscreteKernel(
		vertices,
		edges,
		edgeProps.edgeTangents,
		alpha,
		beta,
		disjointPairs
	);
	
	// Calculate energy based on configuration
	const useSubvertices = get(config).useSubverticesInEnergy;
	const discreteEnergy = useSubvertices
		? calculateDiscreteEnergyWithSubvertices(vertices, edges, alpha, beta, disjointPairs)
		: calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
	
	// Update subvertices based on new vertex positions
	const newSubvertices = generateSubvertices(vertices, edges);
	subvertices.set(newSubvertices);

	return {
		kernelMatrix,
		discreteEnergy,
		edgeProps
	};
}