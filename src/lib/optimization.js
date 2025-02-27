// src/lib/optimization.js
import {
    calculateDifferential,
    calculateEdgeProperties,
    calculateDiscreteEnergy,
    calculateDiscreteEnergyWithSubvertices
} from '$lib/energyCalculations';
import { computePreconditionedGradient, build_A_bar_2D } from '$lib/innerProduct';
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from '$lib/stores';
import { updateKernelState } from '$lib/graphState';

// Import the simplified constraint functionality
import { 
    applyConstraints, 
    projectGradientOntoConstraints
} from '$lib/constraints';

// Configuration for available gradient methods
const GRADIENT_METHODS = {
    L2: 'l2',
    PRECONDITIONED: 'preconditioned'
};

// Default optimization settings
let optimizationConfig = {
    gradientMethod: GRADIENT_METHODS.PRECONDITIONED,
    constraints: {
        maintainEdgeLengths: true,    // Keep individual edge lengths constant
        totalLength: false,           // Control total curve length 
        targetTotalLength: null,      // Target length value (null = use initial)
        barycenter: true,            // Fix curve barycenter
        barycenterTarget: [100, 100]      // Target barycenter position
    }
};

/**
 * Compute gradient based on selected method
 * @param {string} method - Gradient method to use
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {number} alpha - Energy parameter
 * @param {number} beta - Energy parameter
 * @param {Array} disjointPairs - Disjoint edge pairs
 * @returns {Object} - Computed gradient and direction
 */
function computeGradient(method, vertices, edges, alpha, beta, disjointPairs) {
    // Calculate differential (common to all methods)
    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    
    // Based on method, compute the appropriate gradient
    switch(method) {
        case GRADIENT_METHODS.PRECONDITIONED:
            try {
                const { edgeTangents } = calculateEdgeProperties(vertices, edges);
                const gradient = computePreconditionedGradient(
                    vertices,
                    edges,
                    edgeTangents,
                    alpha,
                    beta,
                    differential
                );
                // Descent direction for preconditioned gradient is negative gradient
                const direction = gradient.map(([gx, gy]) => [-gx, -gy]);
                return { gradient, direction, differential };
            } catch (e) {
                console.warn('Preconditioned gradient failed, falling back to L2 gradient:', e);
                // Fall through to L2 method
            }
            
        case GRADIENT_METHODS.L2:
        default:
            // For L2 gradient, the gradient is the differential and direction is negative gradient
            const gradient = differential.map(([dx, dy]) => [dx, dy]);
            const direction = gradient.map(([gx, gy]) => [-gx, -gy]);
            return { gradient, direction, differential };
    }
}

/**
 * Take a fixed-size step in the given direction
 * @param {Array} vertices - Current vertex positions
 * @param {Array} direction - Step direction for each vertex
 * @param {number} stepSize - Step size
 * @returns {Array} - New vertex positions
 */
function takeFixedStep(vertices, direction, stepSize) {
    return vertices.map((vertex, i) => [
        vertex[0] + direction[i][0] * stepSize,
        vertex[1] + direction[i][1] * stepSize
    ]);
}

/**
 * Perform backtracking line search to find optimal step size
 * @param {Array} vertices - Current vertex positions
 * @param {Array} edges - Edge connections
 * @param {Array} direction - Step direction
 * @param {Array} differential - Energy differential
 * @param {number} alpha - Energy parameter
 * @param {number} beta - Energy parameter
 * @param {Array} disjointPairs - Disjoint edge pairs
 * @param {Object} lineSearchSettings - Line search settings
 * @param {Object} constraints - Constraints configuration
 * @param {Object} additionalData - Additional data for constraints
 * @returns {Array} - New vertex positions
 */
function performLineSearch(
    vertices, 
    edges, 
    direction, 
    differential, 
    alpha, 
    beta, 
    disjointPairs, 
    lineSearchSettings, 
    constraints, 
    additionalData
) {
    const directionFlat = direction.flat();
    const differentialFlat = differential.flat();
    const slope = math.dot(differentialFlat, directionFlat);
    
    const currentEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    const { initialStepSize, decay, sufficientDecrease, maxIterations } = lineSearchSettings;
    
    let t = initialStepSize;
    let newVertices;
    
    for (let i = 0; i < maxIterations; i++) {
        // Take step
        const stepVertices = takeFixedStep(vertices, direction, t);
        
        // Apply constraints
        newVertices = applyConstraints(stepVertices, edges, constraints, additionalData);
        
        // Check if step is acceptable
        const newEnergy = calculateDiscreteEnergy(newVertices, edges, alpha, beta, disjointPairs);
        if (newEnergy <= currentEnergy + sufficientDecrease * t * slope) {
            console.log(`Line search converged at t=${t}, energy change: ${currentEnergy - newEnergy}`);
            return newVertices;
        }
        
        // Reduce step size
        t *= decay;
    }
    
    console.warn('Line search did not converge, using smallest step size');
    return newVertices;
}

/**
 * Take a single gradient descent step
 * @param {Array} vertices - Current vertex positions
 * @param {Array} edges - Edge connections
 * @param {number} alpha - Energy parameter
 * @param {number} beta - Energy parameter
 * @param {Array} disjointPairs - Disjoint edge pairs
 * @param {Array} initialEdgeLengths - Initial edge lengths
 * @returns {Array} - New vertex positions
 */
function gradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths) {
    // Get configuration
    const useLineSearch = get(config).useLineSearch;
    const constraints = optimizationConfig.constraints;
    
    // 1. Compute gradient and direction
    const { gradient, direction, differential } = computeGradient(
        optimizationConfig.gradientMethod,
        vertices, 
        edges, 
        alpha, 
        beta, 
        disjointPairs
    );
    
    // 2. Project gradient onto constraint tangent space
    const projectedDirection = projectGradientOntoConstraints(
        direction, 
        vertices, 
        edges, 
        constraints
    );
    
    // 3. Additional data for constraint application
    const additionalData = { 
        initialEdgeLengths,
        alpha,
        beta
    };
    
    // 4. Take step (either fixed or with line search)
    if (useLineSearch) {
        const lineSearchSettings = {
            initialStepSize: get(config).tauInitial,
            decay: get(config).bConst,
            sufficientDecrease: get(config).aConst,
            maxIterations: get(config).maxLineSearch
        };
        
        return performLineSearch(
            vertices, 
            edges, 
            projectedDirection, 
            differential, 
            alpha, 
            beta, 
            disjointPairs,
            lineSearchSettings,
            constraints,
            additionalData
        );
    } else {
        // Fixed step size based on gradient method
        const stepSize = optimizationConfig.gradientMethod === GRADIENT_METHODS.PRECONDITIONED
            ? get(config).precondStepSize
            : get(config).l2StepSize;
            
        const newVertices = takeFixedStep(vertices, projectedDirection, stepSize);
        return applyConstraints(newVertices, edges, constraints, additionalData);
    }
}

/**
 * Create an optimizer for gradient descent
 */
export function createOptimizer(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    maxIterations,
    onUpdate,
    initialEdgeLengths
) {
    if (typeof onUpdate !== 'function') throw new Error('onUpdate must be a function');

    // Initialize target length if using length constraint
    if (optimizationConfig.constraints.totalLength && !optimizationConfig.constraints.targetTotalLength) {
        const totalLength = initialEdgeLengths.reduce((sum, length) => sum + length, 0);
        optimizationConfig.constraints.targetTotalLength = totalLength;
        console.log(`Initialized target total length: ${totalLength}`);
    }

    let currentIteration = 0;
    let intervalId = null;
    let lastEnergy = null;
    let stuckCounter = 0;
    
    const applyPerturbation = get(config).applyPerturbation;
    if (applyPerturbation) {
        lastEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    }

    const optimizer = {
        step: () => {
            if (currentIteration < maxIterations) {
                const newVertices = gradientDescentStep(
                    vertices,
                    edges,
                    alpha,
                    beta,
                    disjointPairs,
                    initialEdgeLengths
                );
                vertices.forEach((v, i) => {
                    v[0] = newVertices[i][0];
                    v[1] = newVertices[i][1];
                });
                
                // Update kernel state including subvertices
                updateKernelState(vertices, edges, alpha, beta, disjointPairs);
                
                if (applyPerturbation) {
                    const newEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
                    const energyChange = newEnergy - lastEnergy;
                    
                    if (Math.abs(energyChange) < (get(config).minEnergyChange || 1e-6)) {
                        stuckCounter++;
                        
                        if (stuckCounter > (get(config).maxStuckIterations || 10)) {
                            console.log('Optimizer stuck, applying random perturbation');
                            applyRandomPerturbation(vertices, get(config).perturbationScale || 0.1);
                            stuckCounter = 0;
                        }
                    } else {
                        stuckCounter = 0;
                    }
                    
                    lastEnergy = newEnergy;
                }
                
                currentIteration++;
                onUpdate();
            } else {
                optimizer.stop();
            }
        },
        start: () => {
            currentIteration = 0;
            stuckCounter = 0;
            
            if (applyPerturbation) {
                lastEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
            }
            
            if (!intervalId) intervalId = setInterval(optimizer.step, 20);
        },
        stop: () => {
            if (intervalId) {
                clearInterval(intervalId);
                intervalId = null;
            }
        },
        // Update gradient method
        setGradientMethod: (method) => {
            if (Object.values(GRADIENT_METHODS).includes(method)) {
                optimizationConfig.gradientMethod = method;
            } else {
                console.warn(`Unknown gradient method: ${method}`);
            }
        },
        // Update constraint settings
        setConstraints: (newConstraints) => {
            optimizationConfig.constraints = {
                ...optimizationConfig.constraints,
                ...newConstraints
            };
        },
        // Get current configuration
        getConfig: () => ({...optimizationConfig})
    };

    return optimizer;
}

function applyRandomPerturbation(vertices, scale) {
    for (const vertex of vertices) {
        vertex[0] += (Math.random() - 0.5) * scale;
        vertex[1] += (Math.random() - 0.5) * scale;
    }
}

// Export constants for use in other files
export const GradientMethods = GRADIENT_METHODS;