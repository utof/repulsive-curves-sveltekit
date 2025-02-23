// src/lib/stores.js
import { writable, derived } from 'svelte/store';

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

// Derived store for discrete energy
export const discreteEnergy = derived(kernelData, ($kernelData) => $kernelData.discreteEnergy);
