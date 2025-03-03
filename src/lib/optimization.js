// src/lib/optimization.js
import {
    calculateDifferential,
    calculateEdgeProperties,
    calculateDiscreteEnergy,
    calculateDiscreteEnergyWithSubvertices
} from '$lib/energyCalculations';
import { computePreconditionedGradient } from '$lib/innerProduct';
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config, optimizationConfig, initialTotalLength, initialEdgeLengths } from '$lib/stores';
import { updateKernelState } from '$lib/graphState';
import { calculateTotalLength } from '$lib/constraints';

// Import the constraint functionality
import { 
    createConstraintData
} from '$lib/constraintJacobian';

// Import constraint projection operations
import {
    projectGradient,
    projectOntoConstraintSet
} from '$lib/constraintProjection';

// Configuration for available gradient methods
export const GradientMethods = {
    L2: 'l2',
    PRECONDITIONED: 'preconditioned'
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
        case GradientMethods.PRECONDITIONED:
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
            
        case GradientMethods.L2:
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
        
        // Get current optimization config
        const optConfig = get(optimizationConfig);
        
        // Apply constraint projection only if there are active constraints
        if (optConfig.useFullConstraintProjection && constraintData.values && constraintData.values.length > 0) {
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
    
    // Get current optimization configuration
    const optConfig = get(optimizationConfig);
    const configValues = get(config);
    
    // Create constraints object from optimization config
    const constraints = {
        barycenter: optConfig.constraints.barycenter.enabled,
        barycenterTarget: optConfig.constraints.barycenter.target,
        length: optConfig.constraints.length.enabled,
        lengthTarget: optConfig.constraints.length.usePercentage ? undefined : optConfig.constraints.length.absoluteValue,
        lengthPercentage: optConfig.constraints.length.usePercentage ? optConfig.constraints.length.percentage : undefined,
        edgeLength: optConfig.constraints.edgeLength.enabled,
        edgeLengthTargets: optConfig.constraints.edgeLength.preserveInitial ? get(initialEdgeLengths) : optConfig.constraints.edgeLength.targets
    };
    
    // Create constraint data for projection
    const constraintData = createConstraintData(vertices, constraints, edges);
    
    // Check if there are any active constraints
    const hasActiveConstraints = constraintData.values && constraintData.values.length > 0;
    
    if (hasActiveConstraints) {
        console.log(`Active constraints: ${Object.entries(constraints)
            .filter(([k, v]) => v === true)
            .map(([k]) => k)
            .join(', ')}`);
        console.log(`Constraint values: [${constraintData.values.join(', ')}]`);
    } else {
        console.log("No active constraints");
    }
    
    // 1. Compute gradient and direction
    const { gradient, direction, differential } = computeGradient(
        optConfig.gradientMethod,
        vertices, 
        edges, 
        alpha, 
        beta, 
        disjointPairs
    );
    
    // 2. Project gradient onto constraint tangent space (Section 5.3.1 of the paper)
    let projectedDirection;
    
    if (optConfig.useFullConstraintProjection && hasActiveConstraints) {
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
    
    // 3. Take step (either fixed or with line search)
    if (optConfig.useLineSearch) {
        const lineSearchSettings = {
            initialStepSize: configValues.tauInitial,
            decay: configValues.bConst,
            sufficientDecrease: configValues.aConst,
            maxIterations: configValues.maxLineSearch
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
        let stepSize = optConfig.gradientMethod === GradientMethods.PRECONDITIONED
            ? optConfig.precondStepSize
            : optConfig.l2StepSize;
        
        // CRITICAL FIX: Add adaptive step size scaling based on direction magnitude
        if (dirMagnitude * stepSize > MAX_STEP_MAGNITUDE) {
            const oldStepSize = stepSize;
            stepSize = MAX_STEP_MAGNITUDE / dirMagnitude;
            console.log(`Reducing step size for stability: ${oldStepSize} -> ${stepSize}`);
        }
        
        console.log(`Taking fixed step with size ${stepSize}`);
        const newVertices = takeFixedStep(vertices, projectedDirection, stepSize);
        
        // Apply constraint projection only if there are active constraints
        if (optConfig.useFullConstraintProjection && hasActiveConstraints) {
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
    
    // Initialize alpha and beta in the optimization config
    optimizationConfig.update(config => ({
        ...config,
        alpha,
        beta
    }));
    
    // Calculate and store the initial total curve length for percentage-based constraints
    const initTotalLength = calculateTotalLength(vertices, edges);
    initialTotalLength.set(initTotalLength);
    console.log(`Initial total curve length: ${initTotalLength}`);
    const initialEdgeLens = initializeEdgeLengths(vertices, edges);
    console.log(`Stored initial edge lengths: [${initialEdgeLens.map(l => l.toFixed(2)).join(', ')}]`);

    const configValues = get(config);
    const applyPerturbation = configValues.applyPerturbation;
    if (applyPerturbation) {
        lastEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    }

    console.log("Creating optimizer with parameters:");
    console.log(`- Gradient method: ${get(optimizationConfig).gradientMethod}`);
    console.log(`- Use full constraint projection: ${get(optimizationConfig).useFullConstraintProjection}`);
    console.log(`- Edge count: ${edges.length}`);
    console.log("Edges structure:", JSON.stringify(edges));

    const optimizer = {
        step: () => {
            if (currentIteration < maxIterations) {
                console.log(`\n--- Iteration ${currentIteration + 1}/${maxIterations} ---`);
                
                // Update current iteration in optimization config for adaptive strategies
                optimizationConfig.update(config => ({
                    ...config, 
                    currentIteration
                }));
                
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
                    
                    // Get latest alpha and beta from store
                    const { alpha: currentAlpha, beta: currentBeta } = get(optimizationConfig);
                    
                    const newVertices = gradientDescentStep(
                        vertices,
                        edges,
                        currentAlpha,
                        currentBeta,
                        disjointPairs
                    );
                    
                    // Update vertices in place
                    vertices.forEach((v, i) => {
                        v[0] = newVertices[i][0];
                        v[1] = newVertices[i][1];
                    });
                    
                    // Update kernel state including subvertices
                    updateKernelState(vertices, edges, currentAlpha, currentBeta, disjointPairs);
                    
                    if (applyPerturbation) {
                        const newEnergy = calculateDiscreteEnergy(vertices, edges, currentAlpha, currentBeta, disjointPairs);
                        const energyChange = newEnergy - lastEnergy;
                        
                        console.log(`Energy: ${newEnergy.toFixed(6)}, Change: ${energyChange.toFixed(6)}`);
                        
                        const configValues = get(config);
                        if (Math.abs(energyChange) < (configValues.minEnergyChange || 1e-6)) {
                            stuckCounter++;
                            console.log(`Possibly stuck: counter = ${stuckCounter}/${configValues.maxStuckIterations || 10}`);
                            
                            if (stuckCounter > (configValues.maxStuckIterations || 10)) {
                                console.log('Optimizer stuck, applying random perturbation');
                                applyRandomPerturbation(vertices, configValues.perturbationScale || 0.1);
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
                const { alpha: currentAlpha, beta: currentBeta } = get(optimizationConfig);
                lastEnergy = calculateDiscreteEnergy(vertices, edges, currentAlpha, currentBeta, disjointPairs);
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
            if (Object.values(GradientMethods).includes(method)) {
                console.log(`Changing gradient method to ${method}`);
                optimizationConfig.update(config => ({
                    ...config,
                    gradientMethod: method
                }));
            } else {
                console.warn(`Unknown gradient method: ${method}`);
            }
        },
        // Update constraint settings
        setConstraints: (newConstraints) => {
            console.log(`Updating constraints: ${JSON.stringify(newConstraints)}`);
            
            optimizationConfig.update(config => ({
                ...config,
                constraints: {
                    barycenter: {
                        enabled: newConstraints.barycenter ?? config.constraints.barycenter.enabled,
                        target: newConstraints.barycenterTarget ?? config.constraints.barycenter.target
                    },
                    length: {
                        enabled: newConstraints.length ?? config.constraints.length.enabled,
                        usePercentage: newConstraints.lengthPercentage !== undefined,
                        percentage: newConstraints.lengthPercentage ?? config.constraints.length.percentage,
                        absoluteValue: newConstraints.lengthTarget ?? config.constraints.length.absoluteValue
                    },
                    // Add edge length constraints
                    edgeLength: {
                        enabled: newConstraints.edgeLength ?? config.constraints.edgeLength.enabled,
                        preserveInitial: newConstraints.preserveInitialEdgeLengths ?? config.constraints.edgeLength.preserveInitial,
                        targets: newConstraints.edgeLengthTargets ?? config.constraints.edgeLength.targets
                    }
                }
            }));
        },
        // Toggle full constraint projection
        setUseFullConstraintProjection: (useFullProjection) => {
            console.log(`Setting full constraint projection: ${useFullProjection}`);
            optimizationConfig.update(config => ({
                ...config,
                useFullConstraintProjection: !!useFullProjection
            }));
        },
        // Update optimizer settings
        updateSettings: (newSettings) => {
            console.log(`Updating optimizer settings: ${JSON.stringify(newSettings)}`);
            optimizationConfig.update(config => ({
                ...config,
                useLineSearch: newSettings.useLineSearch ?? config.useLineSearch,
                precondStepSize: newSettings.precondStepSize ?? config.precondStepSize,
                l2StepSize: newSettings.l2StepSize ?? config.l2StepSize
            }));
        },
        // Update alpha and beta parameters
        updateAlphaBeta: (newAlpha, newBeta) => {
            console.log(`Updating alpha=${newAlpha}, beta=${newBeta}`);
            optimizationConfig.update(config => ({
                ...config,
                alpha: newAlpha,
                beta: newBeta
            }));
        },
        // Get current configuration
        getConfig: () => get(optimizationConfig),
        
        // Reset iterations counter
        resetIterations: () => {
            currentIteration = 0;
            optimizationConfig.update(config => ({
                ...config,
                currentIteration: 0
            }));
        }
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

export function initializeEdgeLengths(vertices, edges) {
    const edgeLengths = [];
    for (const [v1, v2] of edges) {
        const dx = vertices[v2][0] - vertices[v1][0];
        const dy = vertices[v2][1] - vertices[v1][1];
        edgeLengths.push(Math.sqrt(dx*dx + dy*dy));
    }
    initialEdgeLengths.set(edgeLengths);
    return edgeLengths;
}