// src/lib/constraintJacobian.js
import { get } from 'svelte/store';
import { config } from './stores';
import { 
    calculateBarycenter,
    calculateTotalLength,
    evaluateConstraints, // Add this import
    ConstraintNumerics
} from './constraints';
/**
 * Builds the Jacobian matrix for the barycenter constraint
 * Each constraint component (x and y) has a row in the Jacobian
 *
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {number} scalingFactor - Optional scaling for numerical stability
 * @returns {Array} - Jacobian matrix for barycenter constraint [2 x 2|V|]
 */
export function buildBarycenterJacobian(vertices, edges, scalingFactor) {
    console.log("==== BUILDING BARYCENTER JACOBIAN ====");
    const numVertices = vertices.length;
    
    // For the x-component of barycenter constraint
    const xRow = new Array(numVertices * 2).fill(0);
    
    // For the y-component of barycenter constraint
    const yRow = new Array(numVertices * 2).fill(0);
    
    // Calculate edge lengths
    const edgeLengths = [];
    for (const [v1, v2] of edges) {
        const dx = vertices[v2][0] - vertices[v1][0];
        const dy = vertices[v2][1] - vertices[v1][1];
        const length = Math.sqrt(dx*dx + dy*dy);
        edgeLengths.push(length);
    }
    
    // Calculate total weight for normalization
    const totalWeight = edgeLengths.reduce((sum, len) => {
        if (!isFinite(len) || (ConstraintNumerics.USE_SAFE_DENOMINATORS && len <= ConstraintNumerics.MIN_LENGTH)) {
            return sum;
        }
        return sum + len;
    }, 0);
    
    // Check total weight
    if (totalWeight <= 0 || !isFinite(totalWeight)) {
        console.warn(`Invalid total weight in Jacobian calculation: ${totalWeight}`);
        // Fallback with zeros
        return [xRow, yRow];
    }
    
    // Get scaling factor for better numerical conditioning
    if (scalingFactor === undefined) {
        scalingFactor = ConstraintNumerics.USE_STABILITY_SCALING ? 
            get(config).barycenterScaling || ConstraintNumerics.DEFAULT_BARYCENTER_SCALING : 1.0;
    }
    
    console.log(`Using barycenter Jacobian scaling factor: ${scalingFactor}`);
    
    /*
     * For the barycenter constraint: Φ_barycenter(γ) := ∑_{I∈E} ℓ_I (x_I - x₀)
     * Each vertex i affects the barycenter through:
     * 1. Direct contribution to x_I for all edges I that include vertex i
     * 2. Contribution to ℓ_I for all edges I that include vertex i
     * 
     * However, for simplicity and numerical stability, we use an approximation
     * treating the weights as constant, focusing on the first-order effect
     * of vertex positions on the weighted barycenter.
     */
    
    // Explain the approach for clarity
    console.log("Building Jacobian using midpoint-based approach");
    console.log("Each vertex contributes to the barycenter through edges it's part of");
    
    // For each edge, compute how vertices affect the weighted barycenter
    for (let i = 0; i < edges.length; i++) {
        const [v1, v2] = edges[i];
        
        // Skip edges with invalid lengths
        if (!isFinite(edgeLengths[i]) || 
            (ConstraintNumerics.USE_SAFE_DENOMINATORS && edgeLengths[i] <= ConstraintNumerics.MIN_LENGTH)) {
            console.log(`Skipping edge ${i} with invalid length ${edgeLengths[i]}`);
            continue;
        }
        
        // Weight is normalized by total weight for better conditioning
        // Each vertex contributes half to the midpoint
        const weight = edgeLengths[i] / totalWeight / 2 * scalingFactor;
        
        // Vertex v1 affects x-coordinate of barycenter
        xRow[v1 * 2] += weight;
        
        // Vertex v1 affects y-coordinate of barycenter
        yRow[v1 * 2 + 1] += weight;
        
        // Vertex v2 affects x-coordinate of barycenter
        xRow[v2 * 2] += weight;
        
        // Vertex v2 affects y-coordinate of barycenter
        yRow[v2 * 2 + 1] += weight;
    }
    
    // Log max values for diagnostics
    const xRowMax = Math.max(...xRow.map(Math.abs));
    const yRowMax = Math.max(...yRow.map(Math.abs));
    console.log(`Barycenter Jacobian: x-row max value=${xRowMax.toExponential(4)}, y-row max value=${yRowMax.toExponential(4)}`);
    console.log("=====================================");
    
    return [xRow, yRow];
}

/**
 * Builds the Jacobian matrix for the length constraint
 * One row for the total length constraint
 *
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {number} scalingFactor - Optional scaling for numerical stability
 * @returns {Array} - Jacobian matrix for length constraint [1 x 2|V|]
 */
export function buildLengthJacobian(vertices, edges, scalingFactor) {
    console.log("==== BUILDING LENGTH JACOBIAN ====");
    const numVertices = vertices.length;
    
    // One row for the length constraint
    const lengthRow = new Array(numVertices * 2).fill(0);
    
    // Get scaling factor for better numerical conditioning
    if (scalingFactor === undefined) {
        scalingFactor = ConstraintNumerics.USE_STABILITY_SCALING ? 
            get(config).lengthScaling || ConstraintNumerics.DEFAULT_LENGTH_SCALING : 1.0;
    }
    
    console.log(`Using length Jacobian scaling factor: ${scalingFactor}`);
    
    /*
     * For the length constraint: Φ_length(γ) := L^0 - ∑_{I∈E} ℓ_I
     * The Jacobian is the derivative of this constraint with respect to each vertex
     * For each edge (v1, v2), we compute:
     * 
     * d(length)/d(v1.x) = -dx/length
     * d(length)/d(v1.y) = -dy/length
     * d(length)/d(v2.x) = dx/length
     * d(length)/d(v2.y) = dy/length
     */
    
    // Explain the approach
    console.log("Computing derivatives of length with respect to each vertex component");
    
    // For each edge, calculate derivatives of length with respect to vertices
    for (let i = 0; i < edges.length; i++) {
        const [v1, v2] = edges[i];
        const dx = vertices[v2][0] - vertices[v1][0];
        const dy = vertices[v2][1] - vertices[v1][1];
        const length = Math.sqrt(dx*dx + dy*dy);
        
        // Skip edges with invalid lengths
        if (!isFinite(length) || 
            (ConstraintNumerics.USE_SAFE_DENOMINATORS && length <= ConstraintNumerics.MIN_LENGTH)) {
            console.log(`Skipping edge ${i} with invalid length ${length}`);
            continue;
        }
        
        // Calculate the factor dx/length and dy/length
        // Apply scaling factor to improve numerical conditioning
        const factorX = dx * scalingFactor / length;
        const factorY = dy * scalingFactor / length;
        
        // v1 derivatives (negative because moving v1 away from v2 increases length)
        // Negative sign because constraint is L^0 - ∑ℓ_I
        lengthRow[v1 * 2] -= factorX;
        lengthRow[v1 * 2 + 1] -= factorY;
        
        // v2 derivatives
        lengthRow[v2 * 2] += factorX;
        lengthRow[v2 * 2 + 1] += factorY;
    }
    
    // Log max value for diagnostics
    const lengthRowMax = Math.max(...lengthRow.map(Math.abs));
    const nonZeroEntries = lengthRow.filter(v => v !== 0).length;
    console.log(`Length Jacobian: max value=${lengthRowMax.toExponential(4)}, non-zero entries=${nonZeroEntries}/${lengthRow.length}`);
    console.log("=================================");
    
    return [lengthRow];
}


/**
 * Builds the Jacobian matrix for the edge length constraints
 * Each edge contributes one row to the Jacobian
 *
 * @param {Array} vertices - Vertex positions
 * @param {Array} edges - Edge connections
 * @param {number} scalingFactor - Optional scaling for numerical stability
 * @returns {Array} - Jacobian matrix for edge length constraints [|E| x 2|V|]
 */
export function buildEdgeLengthJacobian(vertices, edges, scalingFactor) {
    console.log("==== BUILDING EDGE LENGTH JACOBIAN ====");
    const numVertices = vertices.length;
    const jacobian = [];
    
    // Get scaling factor for better numerical conditioning
    if (scalingFactor === undefined) {
        scalingFactor = ConstraintNumerics.USE_STABILITY_SCALING ? 
            get(config).edgeLengthScaling || ConstraintNumerics.DEFAULT_LENGTH_SCALING : 1.0;
    }
    
    console.log(`Using edge length Jacobian scaling factor: ${scalingFactor}`);
    
    // For each edge, calculate derivatives of its length with respect to vertices
    for (let edgeIdx = 0; edgeIdx < edges.length; edgeIdx++) {
        const [v1, v2] = edges[edgeIdx];
        const dx = vertices[v2][0] - vertices[v1][0];
        const dy = vertices[v2][1] - vertices[v1][1];
        const length = Math.sqrt(dx*dx + dy*dy);
        
        // Skip edges with invalid lengths
        if (!isFinite(length) || 
            (ConstraintNumerics.USE_SAFE_DENOMINATORS && length <= ConstraintNumerics.MIN_LENGTH)) {
            console.log(`Skipping edge ${edgeIdx} with invalid length ${length}`);
            
            // Add a row of zeros as a placeholder
            jacobian.push(new Array(numVertices * 2).fill(0));
            continue;
        }
        
        // Create a row for this edge's constraint
        const row = new Array(numVertices * 2).fill(0);
        
        // Calculate derivatives (same as for total length but for just one edge)
        // d(length)/d(v1.x) = -dx/length
        // d(length)/d(v1.y) = -dy/length
        // d(length)/d(v2.x) = dx/length
        // d(length)/d(v2.y) = dy/length
        
        // Apply scaling factor
        const factorX = dx * scalingFactor / length;
        const factorY = dy * scalingFactor / length;
        
        // v1 derivatives (negative)
        row[v1 * 2] -= factorX;
        row[v1 * 2 + 1] -= factorY;
        
        // v2 derivatives (positive)
        row[v2 * 2] += factorX;
        row[v2 * 2 + 1] += factorY;
        
        jacobian.push(row);
    }
    
    console.log(`Built edge length Jacobian with ${jacobian.length} rows`);
    console.log("=====================================");
    
    return jacobian;
}


/**
 * Builds the complete Jacobian matrix for all constraints
 * Each constraint contributes one or more rows to the matrix
 *
 * @param {Array} vertices - Vertex positions
 * @param {Object} constraints - Constraint configuration
 * @param {Array} edges - Edge connections
 * @returns {Array} - Jacobian matrix C where C[i][j] = ∂Φᵢ/∂γⱼ
 */
export function buildConstraintJacobian(vertices, constraints, edges) {
    console.log("==== BUILDING COMPLETE CONSTRAINT JACOBIAN ====");
    const jacobian = [];
    
    // Barycenter constraint
    if (constraints.barycenter) {
        console.log("Including barycenter constraint Jacobian");
        const barycenterJacobian = buildBarycenterJacobian(vertices, edges);
        jacobian.push(...barycenterJacobian);
        console.log(`Added ${barycenterJacobian.length} rows for barycenter constraint`);
    }
    
    // Length constraint
    if (constraints.length) {
        console.log("Including length constraint Jacobian");
        const lengthJacobian = buildLengthJacobian(vertices, edges);
        jacobian.push(...lengthJacobian);
        console.log(`Added ${lengthJacobian.length} rows for length constraint`);
    }
    
    if (constraints.edgeLength && constraints.edgeLengthTargets) {
        console.log("Including edge length constraint Jacobian");
        const edgeLengthJacobian = buildEdgeLengthJacobian(vertices, edges);
        jacobian.push(...edgeLengthJacobian);
        console.log(`Added ${edgeLengthJacobian.length} rows for edge length constraints`);
    }
    
    console.log(`Complete Jacobian has ${jacobian.length} rows and ${jacobian[0]?.length || 0} columns`);
    console.log("===========================================");
    
    return jacobian;
}


/**
 * Creates a constraint data object containing all information needed for constraint handling
 *
 * @param {Array} vertices - Vertex positions
 * @param {Object} constraints - Constraint configuration
 * @param {Array} edges - Edge connections
 * @returns {Object} - Constraint data including values and Jacobian
 */
export function createConstraintData(vertices, constraints, edges) {
    console.log("==== CREATING CONSTRAINT DATA ====");
    
    
    // Evaluate constraints
    const values = evaluateConstraints(vertices, constraints, edges);
    
    // Build Jacobian
    const jacobian = buildConstraintJacobian(vertices, constraints, edges);
    
    // Return object with constraint data and functions to re-evaluate
    const constraintData = {
        values,
        jacobian,
        evaluate: (v, c) => evaluateConstraints(v, c, edges),
        buildJacobian: (v, c) => buildConstraintJacobian(v, c, edges)
    };
    
    console.log(`Created constraint data with ${values.length} values and ${jacobian.length} Jacobian rows`);
    console.log("================================");
    
    return constraintData;
}