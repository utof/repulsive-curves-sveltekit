// src/lib/stores.js
import { writable, derived } from 'svelte/store';
import { GradientMethods } from './optimization'; // Import gradient method constants

// Numerical parameters and basic configuration
export const config = writable({
    // Numerical stability parameters
    epsilonStability: 1e-7,
    epsilonKernel: 1e-6,
    finiteDiffH: 1e-4,
    constraintTolerance: 1e-2,
    differentialMethod: 'finiteDifference',
    
    // Line search parameters
    tauInitial: 1.0,
    aConst: 0.1,
    bConst: 0.5,
    maxLineSearch: 20,
    maxConstraintIterations: 3,
    
    // Runtime parameters
    applyPerturbation: false,
    minEnergyChange: 1e-6,
    maxStuckIterations: 10,
    perturbationScale: 0.1,
    
    // Subvertex options
    subvertexGap: 100, 
    useSubverticesInEnergy: false,
});

// Centralized optimization configuration
export const optimizationConfig = writable({
    // Energy parameters
    alpha: 3,
    beta: 6,
    linkAlphaBeta: true,  // When true, maintain beta = 2*alpha
    
    // Gradient method
    gradientMethod: GradientMethods.PRECONDITIONED,
    
    // Step size control
    useLineSearch: false,
    precondStepSize: 0.000005,
    l2StepSize: 100000,
    
    // Constraint settings
    constraints: {
        barycenter: {
            enabled: true,
            target: [300, 300]
        },
        length: {
            enabled: false,
            usePercentage: true,
            percentage: 100,
            absoluteValue: 0
        },
        // Future constraint types can be added here
    },
    
    // Constraint implementation
    useFullConstraintProjection: true,
    
    // Constraint stabilization factors
    barycenterStabilization: 0.01,
    lengthStabilization: 0.01,
    
    // Current iteration (used for adaptive strategies)
    currentIteration: 0
});

// Vertex and edge data
export const vertices = writable([]); // Supervertices
export const edges = writable([]);
export const subvertices = writable([]); // Subvertices
export const initialTotalLength = writable(0);

// Kernel data and energy
export const kernelData = writable({
    kernelMatrix: null,
    discreteEnergy: 0,
    edgeProps: { edgeLengths: [], edgeTangents: [], edgeMidpoints: [] },
    disjointPairs: []
});
export const energyChange = writable(0);
export const previousEnergy = writable(0);
export const discreteEnergy = derived(kernelData, ($kernelData) => $kernelData.discreteEnergy);

// Canvas transformation
export const canvasTransform = writable({
    offsetX: 0,
    offsetY: 0,
    zoom: 1.0
});

// History for undo/redo functionality
export const optimizationHistory = writable([]);

// Session data - remembers user preferences
export const userPreferences = writable({
    showAdvancedOptions: false,
    darkMode: false,
});

// Current date and username for logging
export const currentDateTimeUTC = writable("2025-03-02 15:21:29");
export const currentUser = writable("utof");