// src/lib/stores.js
import { writable, derived } from 'svelte/store';

export const config = writable({
    dim: 2, // Default to 2D, user can switch to 3
    epsilonStability: 1e-7,
    epsilonKernel: 1e-6,
    finiteDiffH: 1e-4,
    constraintTolerance: 1e-7,
    tauInitial: 1.0,
    aConst: 0.1,
    bConst: 0.5,
    maxLineSearch: 20,
    differentialMethod: 'finiteDifference',
    precondStepSize: 20,
    l2StepSize: 100000,
    applyPerturbation: false,
    subvertexGap: 50
});

export const vertices = writable([]);
export const edges = writable([]);
export const subvertices = writable([]);
export const kernelData = writable({
    kernelMatrix: null,
    discreteEnergy: 0,
    edgeProps: { edgeLengths: [], edgeTangents: [], edgeMidpoints: [] },
    disjointPairs: []
});
export const energyChange = writable(0);
export const previousEnergy = writable(0);

export const canvasTransform = writable({
    offsetX: 0,
    offsetY: 0,
    zoom: 1.0
});

export const discreteEnergy = derived(kernelData, ($kernelData) => $kernelData.discreteEnergy);