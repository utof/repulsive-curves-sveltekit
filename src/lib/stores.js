// src/lib/stores.js
import { writable, derived } from 'svelte/store';

export const config = writable({
	epsilonStability: 1e-7,
	epsilonKernel: 1e-6,
	finiteDiffH: 1e-4,
	tauInitial: 1.0,
	aConst: 0.1,
	bConst: 0.5,
	constraintTolerance: 1e-2,
	differentialMethod: 'finiteDifference',
    precondStepSize: 0.1,
    l2StepSize: 100000,
    applyPerturbation: true,
    subvertexGap: 50, // Desired gap distance between subvertices (pixels)
    useSubverticesInEnergy: false, // Whether to include subvertices in energy calculations
    useLineSearch: false, // Whether to use line search for step size optimization
	maxLineSearch: 20,
    maxConstraintIterations: 10,
});

export const vertices = writable([]); // Supervertices
export const edges = writable([]);
export const subvertices = writable([]); // New store for subvertices
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
	offsetX: 0,
	offsetY: 0,
	zoom: 1.0
});

export const discreteEnergy = derived(kernelData, ($kernelData) => $kernelData.discreteEnergy);