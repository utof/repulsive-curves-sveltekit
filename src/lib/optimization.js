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

// Import the constraint functionality
import { 
    createConstraintData
} from '$lib/constraints';

// Import constraint projection operations
import {
    projectGradient,
    projectOntoConstraintSet
} from '$lib/constraintProjection';

// Configuration for available gradient methods
const GRADIENT_METHODS = {
    L2: 'l2',
    PRECONDITIONED: 'preconditioned'
};
const annealing = false;
const gradCap = false;

// Default optimization settings
let optimizationConfig = {
    gradientMethod: GRADIENT_METHODS.PRECONDITIONED,
    constraints: {
        barycenter: true,            // Fix curve barycenter
        barycenterTarget: [300, 300]  // Target barycenter position
    },
    // Control whether to use full constraint projection from the paper
    useFullConstraintProjection: true
};

/**
 * Flatten 2D array of vectors into a single 1D array
 * @param {Array} vectors - 2D array of vectors
 * @returns {Array} - Flattened 1D array
 */
function flattenVectors(vectors) {
    return vectors.flatMap(v => v);
}

/**
 * Reshape a flattened array back into a 2D array of vectors
 * @param {Array} flat - Flattened 1D array
 * @param {number} dim - Dimension of each vector (2 for 2D)
 * @returns {Array} - 2D array of vectors
 */
function reshapeToVectors(flat, dim = 2) {
    const result = [];
    for (let i = 0; i < flat.length; i += dim) {
        result.push(flat.slice(i, i + dim));
    }
    return result;
}

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
    console.log(`Computing gradient using ${method} method`);
    console.log("Edges structure in computeGradient:", JSON.stringify(edges));
    
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
    constraints
) {
    console.log("Performing line search with backtracking");
    console.log("Edges structure in performLineSearch:", JSON.stringify(edges));
    
    // Prepare constraint data for projection
    const constraintData = createConstraintData(vertices, constraints, edges);
    
    // Calculate directional derivative
    const directionFlat = flattenVectors(direction);
    const differentialFlat = flattenVectors(differential);
    const slope = math.dot(differentialFlat, directionFlat);
    
    console.log(`Initial directional derivative (slope): ${slope.toFixed(6)}`);
    
    const currentEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    console.log(`Current energy: ${currentEnergy.toFixed(6)}`);
    
    const { initialStepSize, decay, sufficientDecrease, maxIterations } = lineSearchSettings;
    
    let t = initialStepSize;
    let newVertices = null;
    
    for (let i = 0; i < maxIterations; i++) {
        // Take step
        const stepVertices = takeFixedStep(vertices, direction, t);
        console.log(`Line search iteration ${i+1}/${maxIterations}, step size t=${t.toFixed(6)}`);
        
        // Apply constraint projection only if there are active constraints
        if (optimizationConfig.useFullConstraintProjection && constraintData.values && constraintData.values.length > 0) {
            // Use full constraint projection from the paper (Section 5.3.2)
            console.log("Using full constraint projection in line search");
            
            // Define parameters for projection
            const sobolevParams = { alpha, beta };
            
            try {
                newVertices = projectOntoConstraintSet(
                    stepVertices,
                    edges,  // Key fix: edges passed directly here
                    constraints,
                    constraintData,
                    sobolevParams
                );
            } catch (error) {
                console.error("Full constraint projection failed in line search:", error);
                throw error;
            }
        } else {
            // No active constraints, just use the vertices after step
            console.log("No active constraints, skipping projection");
            newVertices = stepVertices;
        }
        
        // Check if step is acceptable
        const newEnergy = calculateDiscreteEnergy(newVertices, edges, alpha, beta, disjointPairs);
        console.log(`New energy: ${newEnergy.toFixed(6)}, change: ${(currentEnergy - newEnergy).toFixed(6)}`);
        
        if (newEnergy <= currentEnergy + sufficientDecrease * t * slope) {
            console.log(`Line search converged at t=${t}, energy change: ${(currentEnergy - newEnergy).toFixed(6)}`);
            return newVertices;
        }
        
        // Reduce step size
        t *= decay;
    }
    
    console.warn('Line search did not converge, using smallest step size');
    return newVertices || vertices.map(v => [...v]);
}

/**
 * Take a single gradient descent step
 * @param {Array} vertices - Current vertex positions
 * @param {Array} edges - Edge connections
 * @param {number} alpha - Energy parameter
 * @param {number} beta - Energy parameter
 * @param {Array} disjointPairs - Disjoint edge pairs
 * @returns {Array} - New vertex positions
 */
function gradientDescentStep(vertices, edges, alpha, beta, disjointPairs) {
    console.log("Taking gradient descent step");
    console.log("Edges in gradientDescentStep:", edges);
    console.log("Edges structure in gradientDescentStep:", JSON.stringify(edges));
    console.log("Edge count:", edges.length);
    
    // Get configuration
    const useLineSearch = get(config).useLineSearch;
    const constraints = optimizationConfig.constraints;
    
    // Create constraint data for projection
    const constraintData = createConstraintData(vertices, constraints, edges);
    
    // Check if there are any active constraints
    const hasActiveConstraints = constraintData.values && constraintData.values.length > 0;
    
    if (hasActiveConstraints) {
        console.log(`Active constraints: ${Object.keys(constraints).filter(k => constraints[k]).join(', ')}`);
        console.log(`Constraint values: [${constraintData.values.join(', ')}]`);
    } else {
        console.log("No active constraints");
    }
    
    // 1. Compute gradient and direction
    const { gradient, direction, differential } = computeGradient(
        optimizationConfig.gradientMethod,
        vertices, 
        edges, 
        alpha, 
        beta, 
        disjointPairs
    );
    
    // 2. Project gradient onto constraint tangent space (Section 5.3.1 of the paper)
    let projectedDirection;
    
    if (optimizationConfig.useFullConstraintProjection && hasActiveConstraints) {
        // Use full saddle point system from the paper
        console.log("Using full gradient projection");
        const flatDirection = flattenVectors(direction);
        
        // Define parameters for projection
        const sobolevParams = { alpha, beta };
        
        try {
            const projectedFlat = projectGradient(
                flatDirection,
                vertices,
                [...edges],  // Create a fresh copy of edges to avoid any potential mutation
                constraintData,
                sobolevParams
            );
            
            // Reshape back to 2D array
            projectedDirection = reshapeToVectors(projectedFlat);
        } catch (error) {
            console.error("Full gradient projection failed:", error);
            throw new Error("Constraint projection failed: " + error.message);
        }
    } else {
        // No constraints to project against, use original direction
        console.log("No constraints to project against, using original direction");
        projectedDirection = direction;
    }
    
    // Log the magnitude of the projected direction
    const dirMagnitude = math.norm(flattenVectors(projectedDirection));
    console.log(`Projected direction magnitude: ${dirMagnitude.toFixed(6)}`);
    
    // CRITICAL FIX: Limit step size for stability
    const MAX_STEP_MAGNITUDE = 1000.0;
    
    // If direction is too small, we might be at a local minimum
    if (dirMagnitude < 1e-6) {
        console.log("Direction magnitude is very small, possible local minimum reached");
        return vertices.map(v => [...v]);  // Return copy of current vertices
    }
    
    // ADD: Detect oscillation by tracking the angle between consecutive directions
    // This would require storing the previous direction - can be added if needed
    
    // 3. Take step (either fixed or with line search)
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
            constraints
        );
    } else {
        // Fixed step size based on gradient method
        let stepSize = optimizationConfig.gradientMethod === GRADIENT_METHODS.PRECONDITIONED
            ? get(config).precondStepSize
            : get(config).l2StepSize;
        
        // CRITICAL FIX: Add adaptive step size scaling based on direction magnitude
        if (gradCap) {
            if (dirMagnitude * stepSize > MAX_STEP_MAGNITUDE) {
                const oldStepSize = stepSize;
                stepSize = MAX_STEP_MAGNITUDE / dirMagnitude;
                console.log(`Reducing step size for stability: ${oldStepSize} -> ${stepSize}`);
            }
        }
        
        // ANNEALING: Gradually reduce step size over iterations to help convergence
        if (annealing) {
            const currentIteration = optimizationConfig.currentIteration || 0;
            if (currentIteration > 20) {
                const annealingFactor = Math.max(0.2, 1.0 - (currentIteration - 20) / 100);
                stepSize *= annealingFactor;
                console.log(`Annealing step size: ${stepSize} (factor: ${annealingFactor})`);
            }
        }
        
        console.log(`Taking fixed step with size ${stepSize}`);
        const newVertices = takeFixedStep(vertices, projectedDirection, stepSize);
        
        // Apply constraint projection only if there are active constraints
        if (optimizationConfig.useFullConstraintProjection && hasActiveConstraints) {
            console.log("Applying full constraint projection after step");
            
            // Define parameters for projection
            const sobolevParams = { alpha, beta };
            
            try {
                return projectOntoConstraintSet(
                    newVertices,
                    edges,
                    constraints,
                    constraintData,
                    sobolevParams
                );
            } catch (error) {
                console.error("Full constraint projection failed after step:", error);
                throw error;
            }
        } else {
            // No constraints to project against, return vertices as is
            console.log("No constraints to project against, skipping projection");
            return newVertices;
        }
    }
}

/**
 * Create an optimizer for gradient descent
 * @param {Array} vertices - Initial vertex positions
 * @param {Array} edges - Edge connections
 * @param {number} alpha - Energy parameter
 * @param {number} beta - Energy parameter
 * @param {Array} disjointPairs - Disjoint edge pairs
 * @param {number} maxIterations - Maximum iterations
 * @param {Function} onUpdate - Callback after each step
 * @returns {Object} - Optimizer object with methods to control optimization
 */
export function createOptimizer(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    maxIterations,
    onUpdate
) {
    if (typeof onUpdate !== 'function') throw new Error('onUpdate must be a function');

    let currentIteration = 0;
    let intervalId = null;
    let lastEnergy = null;
    let stuckCounter = 0;
    
    const applyPerturbation = get(config).applyPerturbation;
    if (applyPerturbation) {
        lastEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    }

    console.log("Creating optimizer with parameters:");
    console.log(`- Gradient method: ${optimizationConfig.gradientMethod}`);
    console.log(`- Use full constraint projection: ${optimizationConfig.useFullConstraintProjection}`);
    console.log(`- Constraints: ${JSON.stringify(optimizationConfig.constraints)}`);
    console.log(`- Edge count: ${edges.length}`);
    console.log("Edges structure:", JSON.stringify(edges));

    const optimizer = {
        step: () => {
            if (currentIteration < maxIterations) {
                console.log(`\n--- Iteration ${currentIteration + 1}/${maxIterations} ---`);
                
                // Store current iteration for adaptive strategies
                optimizationConfig.currentIteration = currentIteration;
                
                try {
                    // Ensure edges is a proper array before passing it
                    if (!Array.isArray(edges)) {
                        console.error("Edges is not an array:", edges);
                        throw new Error(`Invalid edges (not an array): ${typeof edges}`);
                    }
                    
                    if (edges.length === 0) {
                        console.error("Edges array is empty");
                        throw new Error("Invalid edges (empty array)");
                    }
                    
                    // Verify edge structure
                    edges.forEach((edge, i) => {
                        if (!Array.isArray(edge) || edge.length !== 2) {
                            console.error(`Invalid edge at index ${i}:`, edge);
                            throw new Error(`Invalid edge format at index ${i}: ${JSON.stringify(edge)}`);
                        }
                    });
                    
                    const newVertices = gradientDescentStep(
                        vertices,
                        edges,
                        alpha,
                        beta,
                        disjointPairs
                    );
                    
                    // Update vertices in place
                    vertices.forEach((v, i) => {
                        v[0] = newVertices[i][0];
                        v[1] = newVertices[i][1];
                    });
                    
                    // Update kernel state including subvertices
                    updateKernelState(vertices, edges, alpha, beta, disjointPairs);
                    
                    if (applyPerturbation) {
                        const newEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
                        const energyChange = newEnergy - lastEnergy;
                        
                        console.log(`Energy: ${newEnergy.toFixed(6)}, Change: ${energyChange.toFixed(6)}`);
                        
                        if (Math.abs(energyChange) < (get(config).minEnergyChange || 1e-6)) {
                            stuckCounter++;
                            console.log(`Possibly stuck: counter = ${stuckCounter}/${get(config).maxStuckIterations || 10}`);
                            
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
                } catch (error) {
                    console.error("Error during optimization step:", error);
                    optimizer.stop();
                }
            } else {
                console.log("Maximum iterations reached, stopping optimizer");
                optimizer.stop();
            }
        },
        start: () => {
            console.log("Starting optimization");
            currentIteration = 0;
            stuckCounter = 0;
            
            if (applyPerturbation) {
                lastEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
            }
            
            if (!intervalId) intervalId = setInterval(optimizer.step, 20);
        },
        stop: () => {
            console.log("Stopping optimization");
            if (intervalId) {
                clearInterval(intervalId);
                intervalId = null;
            }
        },
        // Update gradient method
        setGradientMethod: (method) => {
            if (Object.values(GRADIENT_METHODS).includes(method)) {
                console.log(`Changing gradient method to ${method}`);
                optimizationConfig.gradientMethod = method;
            } else {
                console.warn(`Unknown gradient method: ${method}`);
            }
        },
        // Update constraint settings
        setConstraints: (newConstraints) => {
            console.log(`Updating constraints: ${JSON.stringify(newConstraints)}`);
            optimizationConfig.constraints = {
                ...optimizationConfig.constraints,
                ...newConstraints
            };
        },
        // Toggle full constraint projection
        setUseFullConstraintProjection: (useFullProjection) => {
            console.log(`Setting full constraint projection: ${useFullProjection}`);
            optimizationConfig.useFullConstraintProjection = !!useFullProjection;
        },
        // Get current configuration
        getConfig: () => ({...optimizationConfig})
    };

    return optimizer;
}

/**
 * Apply a small random perturbation to vertices
 * @param {Array} vertices - Vertex positions to perturb
 * @param {number} scale - Scale factor for perturbation
 */
function applyRandomPerturbation(vertices, scale) {
    console.log(`Applying random perturbation with scale ${scale}`);
    for (const vertex of vertices) {
        vertex[0] += (Math.random() - 0.5) * scale;
        vertex[1] += (Math.random() - 0.5) * scale;
    }
}

// Export constants for use in other files
export const GradientMethods = GRADIENT_METHODS;