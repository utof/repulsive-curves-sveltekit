// src/lib/stores.js
import { writable, derived } from 'svelte/store';

export const config = writable({
	epsilonStability: 1e-7, // For small distance checks (e.g., prevent division by zero)
	epsilonKernel: 1e-6, // For kernel denominators (prevent singularities)
	finiteDiffH: 1e-4, // Step size for finite differences
	constraintTolerance: 1e-7, // Tolerance for constraint projection
	tauInitial: 1.0, // Initial time step for line search
	aConst: 0.1, // Armijo condition constant (0 < a_const < 0.5)
	bConst: 0.5, // Backtracking reduction factor (0 < b_const < 1)
	maxLineSearch: 20, // Maximum iterations for line search
	differentialMethod: 'finiteDifference' // or 'finiteDifference' or analytical
});

export const vertices = writable([]);
export const edges = writable([]);
export const kernelData = writable({
	kernelMatrix: null,
	discreteEnergy: 0,
	edgeProps: { edgeLengths: [], edgeTangents: [], edgeMidpoints: [] },
	disjointPairs: []
});
export const energyChange = writable(0);
export const previousEnergy = writable(0);

// Store for canvas transformations
export const canvasTransform = writable({
	offsetX: 0, // Initial pan offset X
	offsetY: 0, // Initial pan offset Y
	zoom: 1.0 // Initial zoom level
});

export const discreteEnergy = derived(kernelData, ($kernelData) => $kernelData.discreteEnergy);
