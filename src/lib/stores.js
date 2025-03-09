// src/lib/stores.js
import { writable, derived, get } from 'svelte/store';
import { GradientMethods } from './optimization'; // Import gradient method constants

// Numerical parameters and basic configuration
export const config = writable({
    // Numerical stability parameters
    epsilonStability: 1e-7,
    epsilonKernel: 1e-7,
    finiteDiffH: 1e-7,
    constraintTolerance: 1e-4,
    differentialMethod: 'finiteDifference',
    maxConstraintIterations: 1,

    barycenterScaling: 1,
    lengthScaling: 1,
    
    // Line search parameters
    tauInitial: 1.0,
    aConst: 0.1,
    bConst: 0.5,
    maxLineSearch: 20,
    
    // Runtime parameters
    applyPerturbation: false,
    minEnergyChange: 1e-6,
    maxStuckIterations: 10,
    perturbationScale: 0.1,
    
    // Subvertex options
    subvertexGap: 100, 
    useSubverticesInEnergy: false,

    // 3D settings
    dimension: 3, // Default to 2D, can be set to 3 for 3D mode
    zoomZ: 0.5, // Z-axis scaling factor for visualization
    defaultZ: 0, // Default Z coordinate for new vertices in 3D mode
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
    precondStepSize: 1e+2,
    l2StepSize: 100000,
    
    // Constraint settings
    constraints: {
        barycenter: {
            enabled: false,
            target: get(config).dimension === 3 ? [300, 300, 0] : [300, 300]
        },
        length: {
            enabled: false,
            usePercentage: true,
            percentage: 100,
            absoluteValue: 0
        },
        edgeLength: {
            enabled: false,
            preserveInitial: true,  // When true, use initial edge lengths as targets
            targets: []  // Custom target lengths for each edge (if preserveInitial is false)
        }
        // Future constraint types can be added here
    },
    
    // Constraint implementation
    useFullConstraintProjection: true,
    
    // Constraint stabilization factors DEPRECATED TO-DELETE FROM HERE AND CONTROLS
    barycenterStabilization: 0.01,
    lengthStabilization: 1,
    
    // Current iteration (used for adaptive strategies)
    currentIteration: 0
});

// Vertex and edge data
export const vertices = writable([]); // Supervertices
export const edges = writable([]);
export const subvertices = writable([]); // Subvertices
export const initialTotalLength = writable(0);
export const initialEdgeLengths = writable([]);

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