// src/lib/graphState.js
import {
	calculateEdgeProperties,
	calculateDisjointEdgePairs,
	calculateDiscreteKernel,
	calculateDiscreteEnergy
} from '$lib/energyCalculations';

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
	const discreteEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

	return {
		kernelMatrix,
		discreteEnergy,
		edgeProps,
		disjointPairs
	};
}

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
	const discreteEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

	return {
		kernelMatrix,
		discreteEnergy,
		edgeProps
	};
}
