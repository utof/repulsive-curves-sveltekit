Directory Structure:

└── ./
    ├── src
    │   ├── lib
    │   │   ├── energyCalculations.js
    │   │   ├── graphDrawing.js
    │   │   ├── graphstate.js
    │   │   ├── graphUtils.js
    │   │   ├── innerProduct.js
    │   │   ├── interaction.js
    │   │   ├── optimization.js
    │   │   └── stores.js
    │   └── routes
    │       └── +page.svelte
    └── testdiscrete.js



---
File: /src/lib/energyCalculations.js
---

// src/lib/energyCalculations.js
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from '$lib/stores';

let logging = false;

export function calculateEdgeProperties(vertices, edges) {
	const edgeLengths = [];
	const edgeTangents = [];
	const edgeMidpoints = [];

	for (const edge of edges) {
		const v1 = vertices[edge[0]];
		const v2 = vertices[edge[1]];

		const dx = v2[0] - v1[0];
		const dy = v2[1] - v1[1];
		const length = Math.sqrt(dx * dx + dy * dy);
		edgeLengths.push(length);

		const unitTangent = length > 0 ? [dx / length, dy / length] : [0, 0];
		edgeTangents.push(unitTangent);

		const midpoint = [
			isNaN(v1[0]) || isNaN(v2[0]) ? 0 : (v1[0] + v2[0]) / 2,
			isNaN(v1[1]) || isNaN(v2[1]) ? 0 : (v1[1] + v2[1]) / 2
		];
		edgeMidpoints.push(midpoint);

		if (logging) {
			console.log(
				`Edge [${edge[0]}, ${edge[1]}]: length = ${length}, tangent = ${unitTangent}, midpoint = ${midpoint}`
			);
		}
	}

	if (logging) {
		console.log('Edge lengths:', edgeLengths);
		console.log('Unit tangents:', edgeTangents);
		console.log('Midpoints:', edgeMidpoints);
	}

	return { edgeLengths, edgeTangents, edgeMidpoints };
}

export function tangentPointKernel(p, q, T, alpha, beta) {
	// Ensure inputs are properly converted to matrices
	const p_ = math.matrix(p);
	const q_ = math.matrix(q);
	const T_ = math.matrix(T);
	const epsilon = get(config).epsilonKernel;

	// Calculate the difference vector
	const diff = math.subtract(p_, q_);
	const diffNorm = math.norm(diff) + epsilon; // Prevent division by zero
	const cross2D = T_.get([0]) * diff.get([1]) - T_.get([1]) * diff.get([0]); // 2D cross product (determinant)

	const numerator = Math.pow(Math.abs(cross2D), alpha);
	const denominator = Math.pow(diffNorm, beta);
	const result = numerator / denominator;

	// Check for NaN or Infinity
	if (!isFinite(result)) {
		console.warn('Invalid kernel result:', result, 'from inputs:', p, q, T, 'with cross2D:', cross2D, 'diffNorm:', diffNorm);
		return 0;
	}

	return result;
}

export function calculateDisjointEdgePairs(edges) {
	const numEdges = edges.length;
	const disjointPairs = [];

	for (let i = 0; i < numEdges; i++) {
		disjointPairs[i] = [];
		for (let j = 0; j < numEdges; j++) {
			if (i === j) continue;

			const edge1 = edges[i];
			const edge2 = edges[j];

			if (
				edge1[0] !== edge2[0] &&
				edge1[0] !== edge2[1] &&
				edge1[1] !== edge2[0] &&
				edge1[1] !== edge2[1]
			) {
				disjointPairs[i].push(j);
			}
		}
	}
	console.log('Calculated disjointPairs:', disjointPairs);
	return disjointPairs;
}

export function calculateDiscreteKernel(vertices, edges, edgeTangents, alpha, beta, disjointPairs) {
	const numEdges = edges.length;
	const kernelMatrix = math.zeros(numEdges, numEdges);

	if (!disjointPairs || !Array.isArray(disjointPairs) || disjointPairs.length === 0) {
		console.warn('No disjoint pairs found, returning zero kernel matrix');
		return kernelMatrix;
	}

	for (let i = 0; i < numEdges; i++) {
		if (!disjointPairs[i]) {
			console.warn(`No disjoint pairs for edge ${i}`);
			continue;
		}

		for (const j of disjointPairs[i]) {
			if (i < edges.length && j < edges.length) {
				let sum = 0;
				const combinations = [
					[vertices[edges[i][0]], vertices[edges[j][0]]],
					[vertices[edges[i][0]], vertices[edges[j][1]]],
					[vertices[edges[i][1]], vertices[edges[j][0]]],
					[vertices[edges[i][1]], vertices[edges[j][1]]]
				];

				for (const [p, q] of combinations) {
					sum += tangentPointKernel(p, q, edgeTangents[i], alpha, beta);
				}
				kernelMatrix.set([i, j], sum / 4);
				kernelMatrix.set([j, i], sum / 4); // Keep symmetry!
			} else {
				console.warn(
					'Invalid edge index:',
					i,
					j,
					'disjointPairs:',
					disjointPairs,
					'edges.length',
					edges.length
				);
			}
		}
	}
	return kernelMatrix;
}

export function calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs) {
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
	const kernelMatrix = calculateDiscreteKernel(
		vertices,
		edges,
		edgeTangents,
		alpha,
		beta,
		disjointPairs
	);

	let totalEnergy = 0;
	const numEdges = edges.length;

	for (let i = 0; i < numEdges; i++) {
		for (const j of disjointPairs[i]) {
			if (i < edges.length && j < edges.length) {
				const kernelValue = kernelMatrix.get([i, j]);
					totalEnergy += kernelValue * edgeLengths[i] * edgeLengths[j];
			}
		}
	}
	return totalEnergy / 2; // Divide by 2 because of symmetry
}

function calculateAnalyticalDifferential(vertices, edges, alpha, beta, disjointPairs) {
    const numVertices = vertices.length;
    const differential = Array(numVertices).fill().map(() => [0, 0]);
    const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);

    for (let p = 0; p < numVertices; p++) {
        let deriv_p = [0, 0];
        const adjacentEdges = edges.filter(([v0, v1]) => v0 === p || v1 === p);

        for (const edgeI of adjacentEdges) {
            const I = edges.indexOf(edgeI);
            const i = edgeI[0] === p ? 0 : 1;
            const i1 = edgeI[i];
            const i2 = edgeI[(i + 1) % 2];
            const l_I = edgeLengths[I];
            const T_I = edgeTangents[I];

            for (const J of disjointPairs[I]) {
                const l_J = edgeLengths[J];
                const T_J = edgeTangents[J];

                for (let j = 0; j < 2; j++) {
                    const j1 = edges[J][j];
                    const p_i1 = vertices[i1];
                    const p_i2 = vertices[i2];
                    const p_j1 = vertices[j1];

                    // Terms from loss_derivative.cpp adapted for 2D
                    const cross_term = [
                        (p_i2[0] - p_j1[0]) * T_I[1] - (p_i2[1] - p_j1[1]) * T_I[0],
                        (p_i1[0] - p_j1[0]) * T_I[1] - (p_i1[1] - p_j1[1]) * T_I[0]
                    ];
                    const cross_norm = Math.sqrt(cross_term[0] * cross_term[0] + cross_term[1] * cross_term[1]);
                    const denom_diff_i1_j1 = [p_i1[0] - p_j1[0], p_i1[1] - p_j1[1]];
                    const denom_diff_i2_j1 = [p_i2[0] - p_j1[0], p_i2[1] - p_j1[1]];
                    const denom_norm_i1_j1 = Math.sqrt(denom_diff_i1_j1[0] * denom_diff_i1_j1[0] + denom_diff_i1_j1[1] * denom_diff_i1_j1[1]);
                    const denom_norm_i2_j1 = Math.sqrt(denom_diff_i2_j1[0] * denom_diff_i2_j1[0] + denom_diff_i2_j1[1] * denom_diff_i2_j1[1]);

                    // Analytical derivative terms (simplified for 2D)
                    const term1 = (1 - alpha) * Math.pow(l_I, -alpha - 1) * [p_i1[0] - p_i2[0], p_i1[1] - p_i2[1]] * Math.pow(cross_norm, alpha) * Math.pow(denom_norm_i1_j1, -beta);
                    // Add more terms as needed from loss_derivative.cpp
                    deriv_p[0] += 0.25 * l_J * term1[0];
                    deriv_p[1] += 0.25 * l_J * term1[1];
                }
            }
        }
        differential[p] = deriv_p;
    }
    return differential;
}

export function calculateDifferential(vertices, edges, alpha, beta, disjointPairs) {
    const method = get(config).differentialMethod;
    if (method === 'finiteDifference') {
        return calculateDifferentialFiniteDifference(vertices, edges, alpha, beta, disjointPairs);
    } else if (method === 'analytical') {
        return calculateAnalyticalDifferential(vertices, edges, alpha, beta, disjointPairs);
    } else {
        throw new Error('Unknown method for differential calculation');
    }
}

function calculateDifferentialFiniteDifference(vertices, edges, alpha, beta, disjointPairs) {
    // Existing finite difference implementation
    const h = get(config).finiteDiffH;
    const numVertices = vertices.length;
    const differential = [];

    const E_original = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

    for (let i = 0; i < numVertices; i++) {
        differential[i] = [0, 0];
        for (let dim = 0; dim < 2; dim++) {
            const vertices_perturbed = vertices.map((v) => [...v]);
            vertices_perturbed[i][dim] += h;
            const E_perturbed = calculateDiscreteEnergy(
                vertices_perturbed,
                edges,
                alpha,
                beta,
                disjointPairs
            );
            differential[i][dim] = (E_perturbed - E_original) / h;
        }
    }
    console.log('Computed differential:', differential);
    return differential;
}

function calculateL2Gradient(vertices, edges, alpha, beta, disjointPairs) {
	const h = get(config).finiteDiffH; // Use config finiteDiffH
	const numVertices = vertices.length;
	const gradient = [];

	const originalEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

	for (let i = 0; i < numVertices; i++) {
		gradient[i] = [0, 0];

		for (let j = 0; j < 2; j++) {
			const originalValue = vertices[i][j];
			vertices[i][j] += h;

			const newEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

			gradient[i][j] = (newEnergy - originalEnergy) / h;

			vertices[i][j] = originalValue;
		}
	}

	return gradient;
}



---
File: /src/lib/graphDrawing.js
---

// src/lib/graphDrawing.js
import * as math from 'mathjs';
import { drawArrow } from '$lib/graphUtils';
import { calculateDifferential } from '$lib/energyCalculations';
import { canvasTransform } from '$lib/stores';
import { get } from 'svelte/store';

export function drawGraph(
	ctx,
	width,
	height,
	vertices,
	edges,
	edgeProps,
	kernelMatrix,
	alpha,
	beta,
	disjointPairs
) {
	// Use get() to access store values in Svelte 5
	const { offsetX, offsetY, zoom } = get(canvasTransform);

	// Clear the canvas completely in screen coordinates before applying transformations
	ctx.save();
	ctx.clearRect(0, 0, width, height);

	// Apply transformations to draw in world coordinates
	ctx.scale(zoom, zoom);
	ctx.translate(offsetX / zoom, offsetY / zoom);

	drawEdges(ctx, vertices, edges, kernelMatrix);
	drawVertices(ctx, vertices, edges, alpha, beta, disjointPairs);
	drawMidpoints(ctx, edges, edgeProps);

	// Restore context after drawing
	ctx.restore();
}

function drawEdges(ctx, vertices, edges, kernelMatrix) {
	const { zoom } = get(canvasTransform);

	if (kernelMatrix && math.isMatrix(kernelMatrix) && kernelMatrix.size()[0] === edges.length) {
		const maxKernelValue = math.max(kernelMatrix) || 1; // Avoid division by zero

		edges.forEach((edge, i) => {
			const totalKernel = edges.reduce((sum, _, j) => {
				return i !== j ? sum + kernelMatrix.get([i, j]) : sum;
			}, 0);
			const avgKernel = totalKernel / (edges.length - 1 || 1);
			const normalizedValue = avgKernel / maxKernelValue;

			const blue = Math.round(255 * (1 - normalizedValue));
			const red = Math.round(255 * normalizedValue);

			ctx.strokeStyle = `rgb(${red}, 0, ${blue})`;
			ctx.lineWidth = (1 + normalizedValue * 4) / zoom;

			ctx.beginPath();
			ctx.moveTo(vertices[edge[0]][0], vertices[edge[0]][1]);
			ctx.lineTo(vertices[edge[1]][0], vertices[edge[1]][1]);
			ctx.stroke();
		});
	} else {
		// Fallback: Draw edges with default style
		ctx.strokeStyle = 'black';
		ctx.lineWidth = 1 / zoom;
		edges.forEach((edge) => {
			ctx.beginPath();
			ctx.moveTo(vertices[edge[0]][0], vertices[edge[0]][1]);
			ctx.lineTo(vertices[edge[1]][0], vertices[edge[1]][1]);
			ctx.stroke();
		});
	}
}

// src/lib/graphDrawing.js
function drawVertices(ctx, vertices, edges, alpha, beta, disjointPairs) {
    const { offsetX, offsetY, zoom } = get(canvasTransform);
    const gradient = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);

    vertices.forEach((vertex, i) => {
        // Draw circle in world coordinates
        ctx.beginPath();
        ctx.arc(vertex[0], vertex[1], 5 / zoom, 0, 2 * Math.PI);
        ctx.fillStyle = 'blue';
        ctx.fill();

        // Calculate screen coordinates
        const screenX = vertex[0] * zoom + offsetX;
        const screenY = vertex[1] * zoom + offsetY;

        // Draw text and arrow in screen coordinates
        ctx.save();
        ctx.setTransform(1, 0, 0, 1, 0, 0); // Reset transformation

        // Draw vertex number
        ctx.fillStyle = 'black';
        ctx.font = '12px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        ctx.fillText(i.toString(), screenX, screenY - 10);

        // Draw gradient arrow
        const gradX = -gradient[i][0];
        const gradY = -gradient[i][1];
        const magnitude = Math.sqrt(gradX * gradX + gradY * gradY) || 1e-6;
        const screenGradX = gradX * zoom;
        const screenGradY = gradY * zoom;
        const screenMagnitude = Math.sqrt(screenGradX * screenGradX + screenGradY * screenGradY);
        const arrowLength = 20; // pixels
        const dirX = screenGradX / screenMagnitude;
        const dirY = screenGradY / screenMagnitude;
        ctx.strokeStyle = 'purple';
        ctx.lineWidth = 2;
        drawArrow(ctx, screenX, screenY, dirX, dirY, arrowLength);

        ctx.restore();
    });
}

function drawMidpoints(ctx, edges, edgeProps) {
    const { offsetX, offsetY, zoom } = get(canvasTransform);

    if (!edgeProps || !edgeProps.edgeMidpoints || edgeProps.edgeMidpoints.length !== edges.length) {
        return;
    }

    edges.forEach((edge, i) => {
        const midpoint = edgeProps.edgeMidpoints[i];
        const length = edgeProps.edgeLengths[i];
        const tangent = edgeProps.edgeTangents[i];

        if (!midpoint || !length || !tangent) return;

        // Draw midpoint circle in world coordinates
        ctx.fillStyle = 'red';
        ctx.beginPath();
        ctx.arc(midpoint[0], midpoint[1], 2 / zoom, 0, 2 * Math.PI);
        ctx.fill();

        // Calculate screen coordinates
        const screenX = midpoint[0] * zoom + offsetX;
        const screenY = midpoint[1] * zoom + offsetY;

        // Draw text and arrow in screen coordinates
        ctx.save();
        ctx.setTransform(1, 0, 0, 1, 0, 0); // Reset transformation

        ctx.fillStyle = 'blue';
        ctx.font = '10px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        ctx.fillText(`${i}L ${length.toFixed(6)}`, screenX, screenY - 15);
        ctx.fillText(`T ${tangent.map((t) => t.toFixed(1))}`, screenX, screenY - 5);
        ctx.fillText(`${edge[0]}, ${edge[1]}`, screenX, screenY + 15);

        ctx.strokeStyle = 'green';
        ctx.lineWidth = 1.5;
        const tangentScreenX = tangent[0] * zoom;
        const tangentScreenY = tangent[1] * zoom;
        const tangentMagnitude = Math.sqrt(tangentScreenX * tangentScreenX + tangentScreenY * tangentScreenY) || 1;
        const dirX = tangentScreenX / tangentMagnitude;
        const dirY = tangentScreenY / tangentMagnitude;
        drawArrow(ctx, screenX, screenY, dirX, dirY, 20);

        ctx.restore();
    });
}

export function drawKernelMatrix(kernelCanvas, kernelMatrix) {
	if (!kernelCanvas || !kernelMatrix || !math.isMatrix(kernelMatrix)) return;

	const ctx = kernelCanvas.getContext('2d');
	const size = kernelMatrix.size();
	const padding = 50;
	const maxWidth = Math.min(window.innerWidth / 2, 800);
	const maxContentWidth = maxWidth - padding * 2;

	const cellSize = Math.min(
		50,
		Math.floor(maxContentWidth / size[1]),
		Math.floor(maxContentWidth / size[0])
	);

	const canvasWidth = Math.min(maxWidth, padding * 2 + size[1] * cellSize);
	const canvasHeight = padding * 2 + size[0] * cellSize;

	kernelCanvas.width = canvasWidth;
	kernelCanvas.height = canvasHeight;
	kernelCanvas.style.width = `${canvasWidth}px`;
	kernelCanvas.style.height = `${canvasHeight}px`;

	ctx.clearRect(0, 0, canvasWidth, canvasHeight);

	const maxVal = math.max(kernelMatrix) || 1; // Avoid division by zero

	for (let i = 0; i < size[0]; i++) {
		for (let j = 0; j < size[1]; j++) {
			const value = kernelMatrix.get([i, j]);
			const intensity = value / maxVal;
			const r = Math.round(255 - intensity * 255);
			const g = Math.round(255 - intensity * 255);
			const b = 255;
			const minOpacity = 0.1;
			const opacity = minOpacity + (1 - minOpacity) * intensity;
			ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${opacity})`;
			ctx.fillRect(padding + j * cellSize, padding + i * cellSize, cellSize, cellSize);
		}
	}

	ctx.strokeStyle = 'black';
	ctx.lineWidth = 1;
	for (let i = 0; i <= size[0]; i++) {
		ctx.beginPath();
		ctx.moveTo(padding, padding + i * cellSize);
		ctx.lineTo(padding + size[1] * cellSize, padding + i * cellSize);
		ctx.stroke();
	}
	for (let j = 0; j <= size[1]; j++) {
		ctx.beginPath();
		ctx.moveTo(padding + j * cellSize, padding);
		ctx.lineTo(padding + j * cellSize, padding + size[0] * cellSize);
		ctx.stroke();
	}

	ctx.fillStyle = 'black';
	ctx.font = '12px Arial';
	ctx.textAlign = 'center';
	for (let j = 0; j < size[1]; j++) {
		ctx.fillText(j.toString(), padding + j * cellSize + cellSize / 2, padding - 10);
	}

	ctx.textAlign = 'right';
	for (let i = 0; i < size[0]; i++) {
		ctx.fillText(i.toString(), padding - 5, padding + i * cellSize + cellSize / 2);
	}
}



---
File: /src/lib/graphstate.js
---

// src/lib/graphState.js
import {
	calculateEdgeProperties,
	calculateDisjointEdgePairs,
	calculateDiscreteKernel,
	calculateDiscreteEnergy
} from '$lib/energyCalculations';

export function initializeKernelState(vertices, edges, alpha, beta) {
	const disjointPairs = calculateDisjointEdgePairs(edges);
	const edgeProps = calculateEdgeProperties(vertices, edges);
	const kernelMatrix = calculateDiscreteKernel(
		vertices,
		edges,
		edgeProps.edgeTangents,
		alpha,
		beta,
		disjointPairs
	);
	const discreteEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

	return {
		kernelMatrix,
		discreteEnergy,
		edgeProps,
		disjointPairs
	};
}

export function updateKernelState(vertices, edges, alpha, beta, disjointPairs) {
	const edgeProps = calculateEdgeProperties(vertices, edges);
	const kernelMatrix = calculateDiscreteKernel(
		vertices,
		edges,
		edgeProps.edgeTangents,
		alpha,
		beta,
		disjointPairs
	);
	const discreteEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

	return {
		kernelMatrix,
		discreteEnergy,
		edgeProps
	};
}



---
File: /src/lib/graphUtils.js
---

/**
 * Calculate color based on vertex energy
 * @param {number} vertexEnergy - Energy value for the vertex
 * @param {Array} vertexData - Array of vertex data containing energy values
 * @returns {string} RGB color string
 */
export function getColor(vertexEnergy, vertexData) {
	const maxEnergy = Math.max(...vertexData.map((v) => v.energy), 1e-9); // Avoid division by zero
	const normalizedEnergy = vertexEnergy / maxEnergy;
	const r = Math.floor(255 * normalizedEnergy);
	const g = 0;
	const b = Math.floor(255 * (1 - normalizedEnergy));
	return `rgb(${r}, ${g}, ${b})`;
}

/**
 * Calculate radius based on vertex energy
 * @param {number} vertexEnergy - Energy value for the vertex
 * @param {Array} vertexData - Array of vertex data containing energy values
 * @returns {number} Radius value
 */
export function getRadius(vertexEnergy, vertexData) {
	const maxEnergy = Math.max(...vertexData.map((v) => v.energy), 1e-9); // Avoid division by zero
	const normalizedEnergy = vertexEnergy / maxEnergy;
	return 2 + 8 * normalizedEnergy;
}

export function drawArrow(ctx, fromX, fromY, dirX, dirY, length) {
	const headLength = 7;
	const headAngle = Math.PI / 6;

	const toX = fromX + dirX * length;
	const toY = fromY + dirY * length;

	// Draw arrow line
	ctx.beginPath();
	ctx.moveTo(fromX, fromY);
	ctx.lineTo(toX, toY);
	ctx.stroke();

	// Draw arrow head
	const angle = Math.atan2(dirY, dirX);
	ctx.beginPath();
	ctx.moveTo(toX, toY);
	ctx.lineTo(
		toX - headLength * Math.cos(angle - headAngle),
		toY - headLength * Math.sin(angle - headAngle)
	);
	ctx.moveTo(toX, toY);
	ctx.lineTo(
		toX - headLength * Math.cos(angle + headAngle),
		toY - headLength * Math.sin(angle + headAngle)
	);
	ctx.stroke();
}

export function generateRandomGraph(width, height) {
	const numVertices = Math.floor(Math.random() * 10) + 5;
	const vertices = [];
	const edges = [];

	for (let i = 0; i < numVertices; i++) {
		vertices.push([
			Math.random() * width, // x within canvas width
			Math.random() * height // y within canvas height
		]);
	}

	const edgeSet = new Set();
	for (let i = 0; i < numVertices; i++) {
		for (let j = i + 1; j < numVertices; j++) {
			if (Math.random() < 0.1) {
				// 5% chance of an edge
				const edge = [i, j];
				const edgeString = `${i}-${j}`;
				if (!edgeSet.has(edgeString)) {
					edges.push(edge);
					edgeSet.add(edgeString);
				}
			}
		}
	}
	return { vertices, edges };
}

export function generateBipartiteGraph(width, height) {
	// 3x3 bipartite graph: vertices split into two sets (left and right), connected across
	const vertices = [];
	const edges = [];

	// Left set (3 vertices)
	for (let i = 0; i < 3; i++) {
		vertices.push([width * 0.3, height * (0.2 + i * 0.3)]);
	}

	// Right set (3 vertices)
	for (let i = 0; i < 3; i++) {
		vertices.push([width * 0.7, height * (0.2 + i * 0.3)]);
	}

	// Connect each left vertex to all right vertices (bipartite)
	for (let i = 0; i < 3; i++) {
		for (let j = 0; j < 3; j++) {
			edges.push([i, 3 + j]); // Connect left (0-2) to right (3-5)
		}
	}

	return { vertices, edges };
}

export function generate2x3BipartiteGraph(width, height) {
	// 2x3 bipartite graph: vertices split into two sets (left and right), connected across
	const vertices = [];
	const edges = [];

	// Left set (2 vertices)
	for (let i = 0; i < 2; i++) {
		vertices.push([width * 0.3, height * (0.3 + i * 0.4)]);
	}

	// Right set (3 vertices)
	for (let i = 0; i < 3; i++) {
		vertices.push([width * 0.7, height * (0.2 + i * 0.3)]);
	}

	// Connect each left vertex to all right vertices (bipartite)
	for (let i = 0; i < 2; i++) {
		for (let j = 0; j < 3; j++) {
			edges.push([i, 2 + j]); // Connect left (0-1) to right (2-4)
		}
	}

	return { vertices, edges };
}



---
File: /src/lib/innerProduct.js
---

// src/lib/innerProduct.js
import * as math from 'mathjs';
import {
	calculateEdgeProperties,
	calculateDisjointEdgePairs,
	tangentPointKernel
} from '$lib/energyCalculations';
import { get } from 'svelte/store';
import { config } from '$lib/stores';

/**
 * Builds the weight matrices W and W0 used in energy calculations.
 * @param {number} alpha - Energy parameter.
 * @param {number} beta - Energy parameter.
 * @param {Array<Array<number>>} edges - Array of edges, where each edge is an array of two vertex indices.
 * @param {Array<Array<number>>} disjointPairs - Array of disjoint edge pairs.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edgeTangents - Array of edge tangent vectors.
 * @param {Array<number>} edgeLengths - Array of edge lengths.
 * @returns {{W: math.Matrix, W0: math.Matrix}} - The weight matrices W and W0.
 */
function build_weights(alpha, beta, edges, disjointPairs, vertices, edgeTangents, edgeLengths) {
	const edge_num = edges.length;
	const s = (beta - 1) / alpha; // Fractional order s from the paper
	const sigma = s - 1; // Correct sigma = s - 1 (Section 4.2.2)

	const W = math.zeros(edge_num, edge_num);
	const W0 = math.zeros(edge_num, edge_num);

	for (let I = 0; I < disjointPairs.length; I++) {
		for (const J of disjointPairs[I]) {
			let elt1 = 0;
			let elt2 = 0;

			for (let a = 0; a < 2; a++) {
				for (let b = 0; b < 2; b++) {
					const i = edges[I][a];
					const j = edges[J][b];
					const p = vertices[i];
					const q = vertices[j];
					const diff = math.subtract(p, q);
					const epsilon = get(config).epsilonStability; // Use config epsilon
					let diff_norm = math.norm(diff) + epsilon; // Prevent division by zero

					const term1 = 1 / Math.pow(diff_norm, 2 * sigma + 1);
					elt1 += term1;

					// Constants for B0 calculation
					const alph = 2;
					const bet = 4;

					const cross = math.det([diff, edgeTangents[I]]); // 2D cross product
					const cross_norm = Math.abs(cross); // Use absolute value for 2D

					const k_numerator = Math.pow(cross_norm, alph);
					const k_denominator = Math.pow(diff_norm, bet);
					const k = k_numerator / k_denominator;

					const term2 = k / Math.pow(diff_norm, 2 * sigma + 1);
					elt2 += term2;
				}
			}

			const w_ij_factor = 0.25 * edgeLengths[I] * edgeLengths[J];
			W.set([I, J], w_ij_factor * elt1);
			W0.set([I, J], w_ij_factor * elt2);
		}
	}

	return { W, W0 };
}

/**
 * Calculates the low-order term matrix B0.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edges - Array of edges.
 * @param {math.Matrix} W0 - The weight matrix W0.
 * @returns {math.Matrix} - The low-order term matrix B0.
 */
function calculateLowOrderTerm(vertices, edges, W0) {
	const numVertices = vertices.length;
	const B0 = math.zeros(numVertices, numVertices);
	const disjointEdges = calculateDisjointEdgePairs(edges);

	for (let I = 0; I < edges.length; I++) {
		for (const J of disjointEdges[I]) {
			const w_IJ_0 = W0.get([I, J]);

			for (let a = 0; a < 2; a++) {
				for (let b = 0; b < 2; b++) {
					const i_a = edges[I][a];
					const i_b = edges[I][b];
					const j_a = edges[J][a];
					const j_b = edges[J][b];

					B0.set([i_a, i_b], B0.get([i_a, i_b]) + 0.25 * w_IJ_0);
					B0.set([j_a, j_b], B0.get([j_a, j_b]) + 0.25 * w_IJ_0);
					B0.set([i_a, j_b], B0.get([i_a, j_b]) - 0.25 * w_IJ_0);
					B0.set([j_a, i_b], B0.get([j_a, i_b]) - 0.25 * w_IJ_0);
				}
			}
		}
	}
	console.log('Low order term B0:', B0.toArray());
	return B0;
}

/**
 * Calculates the high-order term matrix B.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edges - Array of edges.
 * @param {math.Matrix} W - The weight matrix W.
 * @param {Array<number>} edgeLengths - Array of edge lengths.
 * @param {Array<Array<number>>} edgeTangents - Array of edge tangent vectors.
 * @returns {math.Matrix} - The high-order term matrix B.
 */
function calculateHighOrderTerm(vertices, edges, W, edgeLengths, edgeTangents) {
	const numVertices = vertices.length;
	const B = math.zeros(numVertices, numVertices);
	const disjointEdges = calculateDisjointEdgePairs(edges);

	for (let I = 0; I < edges.length; I++) {
		for (const J of disjointEdges[I]) {
			const l_I = edgeLengths[I];
			const l_J = edgeLengths[J];
			const T_I = edgeTangents[I];
			const T_J = edgeTangents[J];
			const w_IJ = W.get([I, J]);
			const dot_TI_TJ = math.dot(T_I, T_J); // 2D dot product

			for (let a = 0; a < 2; a++) {
				for (let b = 0; b < 2; b++) {
					const sign = Math.pow(-1, a + b);
					const i_a = edges[I][a];
					const i_b = edges[I][b];
					const j_a = edges[J][a];
					const j_b = edges[J][b];

					const val_1 = (sign * w_IJ) / (l_I * l_I);
					const val_2 = (sign * w_IJ) / (l_J * l_J);
					const val_3 = (sign * w_IJ * dot_TI_TJ) / (l_I * l_J);

					B.set([i_a, i_b], B.get([i_a, i_b]) + val_1);
					B.set([j_a, j_b], B.get([j_a, j_b]) + val_2);
					B.set([i_a, j_b], B.get([i_a, j_b]) - val_3);
					B.set([j_a, i_b], B.get([j_a, i_b]) - val_3);
				}
			}
		}
	}

	console.log('High order term B:', B.toArray());
	return B;
}

/**
 * Calculates the discrete inner product matrix A and its components.
 * @param {Array<Array<number>>} vertices - Array of vertex coordinates.
 * @param {Array<Array<number>>} edges - Array of edges.
 * @param {number} alpha - Energy parameter.
 * @param {number} beta - Energy parameter.
 * @returns {{A_reg: math.Matrix, B0: math.Matrix, B: math.Matrix}} - The regularized inner product matrix A_reg and component matrices B0 and B.
 */
export function calculateDiscreteInnerProduct(vertices, edges, alpha, beta) {
	console.log('Calculating discrete inner product with alpha:', alpha, 'beta:', beta);
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);

	const { W, W0 } = build_weights(
		alpha,
		beta,
		edges,
		calculateDisjointEdgePairs(edges),
		vertices,
		edgeTangents,
		edgeLengths
	);
	const B = calculateHighOrderTerm(vertices, edges, W, edgeLengths, edgeTangents);
	const B0 = calculateLowOrderTerm(vertices, edges, W0);

	let A = math.add(B0, B);
	console.log('Initial A Matrix (B0 + B):', A.toArray());
	const A_reg = A;

	return { A_reg, B0, B };
}

export function build_A_bar_2D(alpha, beta, vertices, edges) {
	const { edgeLengths, edgeTangents } = calculateEdgeProperties(vertices, edges);
	const disjointEdges = calculateDisjointEdgePairs(edges);
	const numVertices = vertices.length;

	const { W, W0 } = build_weights(
		alpha,
		beta,
		edges,
		disjointEdges,
		vertices,
		edgeTangents,
		edgeLengths
	);
	const B = calculateHighOrderTerm(vertices, edges, W, edgeLengths, edgeTangents);
	const B0 = calculateLowOrderTerm(vertices, edges, W0);
	let A = math.add(B, B0);

	// Regularization to ensure invertibility (use config epsilon)
	const epsilon = get(config).epsilonStability;
	const reg = math.multiply(epsilon, math.identity(numVertices));
	A = math.add(A, reg);

	// Build 2D block-diagonal A_bar
	const A_bar = math.zeros(2 * numVertices, 2 * numVertices);
	A_bar.subset(math.index(math.range(0, numVertices), math.range(0, numVertices)), A);
	A_bar.subset(
		math.index(math.range(numVertices, 2 * numVertices), math.range(numVertices, 2 * numVertices)),
		A
	);

	return { A_bar, B, B0 };
}

export function computePreconditionedGradient(
	vertices,
	edges,
	edgeTangents,
	alpha,
	beta,
	differential
) {
	console.log('Computing preconditioned gradient with differential:', differential);
	const numVertices = vertices.length;
	const { A_bar } = build_A_bar_2D(alpha, beta, vertices, edges);
	const differentialFlat = differential.flat();

	let gradFlat;
	try {
		gradFlat = math.lusolve(A_bar, differentialFlat);
		console.log('Preconditioned gradient (flat):', gradFlat.toArray());
	} catch (e) {
		console.error('Linear solve failed:', e);
		throw new Error('Failed to compute preconditioned gradient due to singular matrix');
	}

	const gradArray = gradFlat.toArray();
	const grad = [];
	for (let i = 0; i < numVertices; i++) {
		grad[i] = [gradArray[i * 2], gradArray[i * 2 + 1]];
	}
	console.log('Gradient:', grad);
	return grad;
}



---
File: /src/lib/interaction.js
---

// src/lib/interaction.js
import { get, writable } from 'svelte/store';
import { canvasTransform } from '$lib/stores';

export function setupInteractions(canvas, vertices, updateFn) {
	let draggingVertex = null;
	let dragOffsetX = 0;
	let dragOffsetY = 0;
	let isDraggingCanvas = false; // For panning the canvas
	let lastMouseX = 0;
	let lastMouseY = 0;
	let isSpacePressed = false; // Track Space key state
	let isCtrlPressed = false; // Track Ctrl key state

	const MIN_ZOOM = 0.1; // Minimum zoom level
	const MAX_ZOOM = 5.0; // Maximum zoom level

	function getWorldCoords(screenX, screenY) {
		const { offsetX, offsetY, zoom } = get(canvasTransform);
		return {
			worldX: (screenX - offsetX) / zoom,
			worldY: (screenY - offsetY) / zoom
		};
	}

	function getScreenCoords(worldX, worldY) {
		const { offsetX, offsetY, zoom } = get(canvasTransform);
		return {
			screenX: worldX * zoom + offsetX,
			screenY: worldY * zoom + offsetY
		};
	}

	function handleMouseDown(event) {
		const rect = canvas.getBoundingClientRect();
		const mouseX = event.clientX - rect.left;
		const mouseY = event.clientY - rect.top;

		// Check if Space is pressed for canvas dragging (panning)
		if (isSpacePressed) {
			isDraggingCanvas = true;
			lastMouseX = mouseX;
			lastMouseY = mouseY;
			return; // Skip vertex dragging if panning
		}

		// Check for vertex dragging (existing behavior)
		const { zoom } = get(canvasTransform);
		for (let i = 0; i < vertices.length; i++) {
			const [vx, vy] = vertices[i];
			const { screenX, screenY } = getScreenCoords(vx, vy);
			const distance = Math.sqrt((mouseX - screenX) ** 2 + (mouseY - screenY) ** 2);
			if (distance <= 10 / zoom) {
				// Increased click radius for better hit detection
				draggingVertex = i;
				const worldCoords = getWorldCoords(mouseX, mouseY);
				dragOffsetX = vx - worldCoords.worldX;
				dragOffsetY = vy - worldCoords.worldY;
				break;
			}
		}
	}

	function handleMouseMove(event) {
		const rect = canvas.getBoundingClientRect();
		const mouseX = event.clientX - rect.left;
		const mouseY = event.clientY - rect.top;

		if (isDraggingCanvas && isSpacePressed) {
			// Panning: Update canvas offset with throttling for performance
			const dx = mouseX - lastMouseX;
			const dy = mouseY - lastMouseY;
			canvasTransform.update((transform) => ({
				...transform,
				offsetX: transform.offsetX + dx,
				offsetY: transform.offsetY + dy
			}));
			lastMouseX = mouseX;
			lastMouseY = mouseY;
			requestAnimationFrame(updateFn); // Use requestAnimationFrame for smoother updates
			return; // Skip vertex dragging if panning
		}

		if (draggingVertex !== null) {
			// Vertex dragging (optimized for performance)
			const worldCoords = getWorldCoords(mouseX, mouseY);
			let newX = worldCoords.worldX + dragOffsetX;
			let newY = worldCoords.worldY + dragOffsetY;

			// Keep vertices within canvas bounds (optional, adjust as needed)
			const { zoom } = get(canvasTransform);
			newX = Math.max(0, Math.min(canvas.width / zoom, newX));
			newY = Math.max(0, Math.min(canvas.height / zoom, newY));

			vertices[draggingVertex] = [newX, newY];
			requestAnimationFrame(updateFn); // Use requestAnimationFrame for smoother updates
		}
	}

	function handleMouseUp() {
		draggingVertex = null;
		isDraggingCanvas = false;
	}

	function handleMouseWheel(event) {
		if (isCtrlPressed) {
			event.preventDefault();
			const rect = canvas.getBoundingClientRect();
			const mouseX = event.clientX - rect.left;
			const mouseY = event.clientY - rect.top;

			// Calculate zoom factor
			const zoomFactor = event.deltaY < 0 ? 1.1 : 0.9; // Zoom in/out by 10%
			canvasTransform.update((transform) => {
				const newZoom = Math.min(MAX_ZOOM, Math.max(MIN_ZOOM, transform.zoom * zoomFactor));
				const worldCoords = getWorldCoords(mouseX, mouseY);
				const newOffsetX = mouseX - worldCoords.worldX * newZoom;
				const newOffsetY = mouseY - worldCoords.worldY * newZoom;
				return {
					offsetX: newOffsetX,
					offsetY: newOffsetY,
					zoom: newZoom
				};
			});

			requestAnimationFrame(updateFn); // Use requestAnimationFrame for smoother updates
		}
	}

	function handleKeyDown(event) {
		if (event.code === 'Space') {
			isSpacePressed = true;
		} else if (event.code === 'ControlLeft' || event.code === 'ControlRight') {
			isCtrlPressed = true;
		}
	}

	function handleKeyUp(event) {
		if (event.code === 'Space') {
			isSpacePressed = false;
			isDraggingCanvas = false; // Stop panning when Space is released
		} else if (event.code === 'ControlLeft' || event.code === 'ControlRight') {
			isCtrlPressed = false;
		}
	}

	// Add event listeners
	canvas.addEventListener('mousedown', handleMouseDown);
	canvas.addEventListener('mousemove', handleMouseMove);
	canvas.addEventListener('mouseup', handleMouseUp);
	canvas.addEventListener('mouseleave', handleMouseUp);
	canvas.addEventListener('wheel', handleMouseWheel, { passive: false });
	document.addEventListener('keydown', handleKeyDown);
	document.addEventListener('keyup', handleKeyUp);

	// Return a cleanup function
	return () => {
		canvas.removeEventListener('mousedown', handleMouseDown);
		canvas.removeEventListener('mousemove', handleMouseMove);
		canvas.removeEventListener('mouseup', handleMouseUp);
		canvas.removeEventListener('mouseleave', handleMouseUp);
		canvas.removeEventListener('wheel', handleMouseWheel);
		document.removeEventListener('keydown', handleKeyDown);
		document.removeEventListener('keyup', handleKeyUp);
	};
}



---
File: /src/lib/optimization.js
---

// src/lib/optimization.js
import {
    calculateDifferential,
    calculateEdgeProperties,
    calculateDiscreteEnergy
} from '$lib/energyCalculations';
import { computePreconditionedGradient } from '$lib/innerProduct';
import * as math from 'mathjs';
import { get } from 'svelte/store';
import { config } from '$lib/stores';

// Configuration toggles
const usePreconditioned = false; // Use preconditioned gradient by default
const applyProjectConstraints = false; // Enable length constraints
const applyBarycenter = true; // Enable barycenter constraint

function l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths) {
    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    const gradient = differential.map(([dx, dy]) => [dx, dy]);
    
    // Normalize the gradient to avoid too large steps
    // const gradNorm = Math.sqrt(gradient.flat().reduce((sum, val) => sum + val * val, 0)) || 1;
    const stepSize = get(config).l2StepSize; // Get step size from config
    
    // Create new vertices by taking a step in the negative gradient direction
    const newVertices = vertices.map((vertex, i) => [
        vertex[0] - stepSize * gradient[i][0] ,
        vertex[1] - stepSize * gradient[i][1]
    ]);

    // Apply constraints if enabled
    let constrainedVertices = [...newVertices];
    if (applyProjectConstraints) {
        constrainedVertices = projectConstraints(constrainedVertices, edges, initialEdgeLengths);
    }
    
    
    if (applyBarycenter) {
        enforceBarycenter(constrainedVertices);
    }

    return constrainedVertices;
}

function preconditionedGradientDescentStep(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    initialEdgeLengths
) {
    const { edgeTangents } = calculateEdgeProperties(vertices, edges);
    
    const differential = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);
    
    // Compute the preconditioned gradient
    let gradient;
    try {
        gradient = computePreconditionedGradient(
            vertices,
            edges,
            edgeTangents,
            alpha,
            beta,
            differential
        );
    } catch (e) {
        console.warn('Preconditioned gradient failed, falling back to L2 gradient:', e);
        return l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths);
    }




    const d = gradient.map(([gx, gy]) => [-gx, -gy]);
    const d_norm = Math.sqrt(d.flat().reduce((sum, val) => sum + val * val, 0)) || 1;
    const d_normalized = d.map(([dx, dy]) => [dx / d_norm, dy / d_norm]);

    const differentialFlat = differential.flat();
    const dFlat = d_normalized.flat();
    const slope = math.dot(differentialFlat, dFlat);
    
    const E_old = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
    
    // Line search parameters
    const a_const = get(config).aConst;
    const b_const = get(config).bConst;
    const max_line_search = get(config).maxLineSearch;
    let t = get(config).tauInitial;
    const stepSize = get(config).precondStepSize; // Get step size from config

    // Perform line search to find a good step size
    let vertices_new = [...vertices];
    for (let i = 0; i < max_line_search; i++) {
        vertices_new = vertices.map((vertex, idx) => [
            vertex[0] + t * d_normalized[idx][0] * stepSize,
            vertex[1] + t * d_normalized[idx][1] * stepSize
        ]);
        
        // Apply constraints
        if (applyProjectConstraints) {
            vertices_new = projectConstraints(vertices_new, edges, initialEdgeLengths);
        }
        
        if (applyBarycenter) {
            enforceBarycenter(vertices_new);
        }
        
        const E_new = calculateDiscreteEnergy(vertices_new, edges, alpha, beta, disjointPairs);
        console.log(`Line search iteration ${i}: t=${t}, E_new=${E_new}, E_old=${E_old}, condition=${E_old + a_const * t * slope}`);
        if (E_new <= E_old + a_const * t * slope) {
            console.log('Line search converged at t=', t, 'with energy reduction:', E_old - E_new);
            return vertices_new;
        }
        
        // Backtracking: reduce step size
        t *= b_const;
    }

    console.warn('Line search did not converge, using smallest step size');
    return vertices_new;
}

// Improved constraint projection to maintain edge lengths
function projectConstraints(
    vertices,
    edges,
    initialEdgeLengths,
    maxIterations = get(config).maxConstraintIterations,
    tolerance = get(config).constraintTolerance
) {
    // Create a deep copy to avoid modifying the input
    const projectedVertices = vertices.map(v => [...v]);
    
    for (let iter = 0; iter < maxIterations; iter++) {
        let maxError = 0;
        
        for (let i = 0; i < edges.length; i++) {
            const [v1Idx, v2Idx] = edges[i];
            const v1 = projectedVertices[v1Idx];
            const v2 = projectedVertices[v2Idx];
            
            // Calculate current edge vector and length
            const dx = v2[0] - v1[0];
            const dy = v2[1] - v1[1];
            const currentLength = Math.sqrt(dx * dx + dy * dy) + get(config).epsilonStability;

            // Calculate error relative to target length
            const targetLength = initialEdgeLengths[i];
            const error = (currentLength - targetLength) / targetLength;
            
            // Update the maximum error
            maxError = Math.max(maxError, Math.abs(error));
            
            // If error is significant, adjust vertex positions
            if (Math.abs(error) > tolerance) {
                const correction = error / 2; // Split the correction between vertices
                const scaleFactor = 1 - correction;
                
                // Calculate midpoint
                const midX = (v1[0] + v2[0]) / 2;
                const midY = (v1[1] + v2[1]) / 2;
                
                // Adjust vertex positions to preserve midpoint
                const halfDx = dx / 2;
                const halfDy = dy / 2;
                
                projectedVertices[v1Idx][0] = midX - halfDx * scaleFactor;
                projectedVertices[v1Idx][1] = midY - halfDy * scaleFactor;
                projectedVertices[v2Idx][0] = midX + halfDx * scaleFactor;
                projectedVertices[v2Idx][1] = midY + halfDy * scaleFactor;
            }
        }
        
        // If all constraints are satisfied to tolerance, exit early
        if (maxError < tolerance) {
            console.log(`Constraint projection converged after ${iter + 1} iterations`);
            break;
        }
    }
    
    return projectedVertices;
}

// Improved barycenter constraint
// Improved barycenter constraint with options
function enforceBarycenter(vertices, options = {}) {
    // Default options
    const {
        targetBarycenter = null,  // If null, maintain original position
        centerAtOrigin = false    // If true, center at (0,0) (original behavior)
    } = options;
    
    // Calculate the current center of mass
    const currentBarycenter = [0, 0];
    const n = vertices.length;
    
    for (const vertex of vertices) {
        currentBarycenter[0] += vertex[0] / n;
        currentBarycenter[1] += vertex[1] / n;
    }
    
    // Determine target position based on options
    let target;
    if (centerAtOrigin) {
        target = [0, 0]; // Original behavior (centers at origin)
    } else if (targetBarycenter) {
        target = targetBarycenter; // Use specified target
    } else {
        target = currentBarycenter; // Maintain current position (no movement)
    }
    
    // Calculate the translation needed
    const dx = target[0] - currentBarycenter[0];
    const dy = target[1] - currentBarycenter[1];
    
    // Apply the translation to all vertices
    for (const vertex of vertices) {
        vertex[0] += dx;
        vertex[1] += dy;
    }
    
    return vertices;
}
export function gradientDescentStep(
    vertices,
    edges,
    alpha,
    beta,
    disjointPairs,
    initialEdgeLengths
) {
    return usePreconditioned
        ? preconditionedGradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths)
        : l2GradientDescentStep(vertices, edges, alpha, beta, disjointPairs, initialEdgeLengths);
}

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

    let currentIteration = 0;
    let intervalId = null;
    let lastEnergy = null;
    let stuckCounter = 0;
    
    // Only initialize energy tracking if perturbation is enabled
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
                
                // Only perform energy calculations and stuck detection if perturbation is enabled
                if (applyPerturbation) {
                    // Calculate new energy
                    const newEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
                    const energyChange = newEnergy - lastEnergy;
                    
                    // Check if we're stuck (minimal energy change)
                    if (Math.abs(energyChange) < get(config).minEnergyChange) {
                        stuckCounter++;
                        
                        // Apply perturbation if stuck for too many iterations
                        if (stuckCounter > get(config).maxStuckIterations) {
                            console.log('Optimizer stuck, applying random perturbation');
                            applyRandomPerturbation(vertices, get(config).perturbationScale);
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
            
            // Initialize energy tracking on start if perturbation is enabled
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
        }
    };

    return optimizer;
}

// Helper function to apply a small random perturbation when optimization gets stuck
function applyRandomPerturbation(vertices, scale) {
    for (const vertex of vertices) {
        vertex[0] += (Math.random() - 0.5) * scale;
        vertex[1] += (Math.random() - 0.5) * scale;
    }
}


---
File: /src/lib/stores.js
---

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
	differentialMethod: 'finiteDifference', // or 'finiteDifference' or analytical
    precondStepSize: 20,
    l2StepSize: 100000,
    applyPerturbation: false
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



---
File: /src/routes/+page.svelte
---

<script>
	import Controls from '$lib/Controls.svelte';
	import { onMount, onDestroy } from 'svelte';
	import { drawGraph, drawKernelMatrix } from '$lib/graphDrawing';
	import { createOptimizer } from '$lib/optimization';
	import {
		generateRandomGraph,
		generateBipartiteGraph,
		generate2x3BipartiteGraph
	} from '$lib/graphUtils';
	import { setupInteractions } from '$lib/interaction';
	import { initializeKernelState, updateKernelState } from '$lib/graphState';
	import {
		vertices,
		edges,
		kernelData,
		energyChange,
		previousEnergy,
		discreteEnergy,
		config,
		canvasTransform
	} from '$lib/stores';
	import { get } from 'svelte/store';

	let graphCanvas;
	let kernelCanvas;
	let graphCtx;
	let optimizer;
	let cleanupInteractions = () => {};
	let isOptimizing = false;
	let graphType = 'bipartite'; // Default to bipartite graph
	// let graphType = 'random'; // Changed to 'random' graph type
	const width = 700;
	const height = 700;
	let alpha = 3; // Default values from paper recommendation
	let beta = 6;
	const maxIterations = 1000;
	let initialEdgeLengths = [];

	onMount(() => {
		graphCtx = graphCanvas.getContext('2d');
		regenerateGraph();
	});

	onDestroy(() => {
		if (optimizer) optimizer.stop();
		cleanupInteractions();
	});

	function updateVisualization() {
		// Ensure a clean redraw by triggering drawGraph
		const updatedKernel = updateKernelState(
			$vertices,
			$edges,
			alpha,
			beta,
			$kernelData.disjointPairs
		);
		$kernelData = { ...updatedKernel, disjointPairs: $kernelData.disjointPairs };
		$energyChange = $discreteEnergy - $previousEnergy;
		$previousEnergy = $discreteEnergy;

		drawGraph(
			graphCtx,
			width,
			height,
			$vertices,
			$edges,
			updatedKernel.edgeProps,
			updatedKernel.kernelMatrix,
			alpha,
			beta,
			$kernelData.disjointPairs
		);

		// Draw kernel matrix (no transformation needed)
		drawKernelMatrix(kernelCanvas, updatedKernel.kernelMatrix);
	}

	function regenerateGraph() {
		$previousEnergy = 0;
		$energyChange = 0;
		let newVertices, newEdges;
		if (graphType === 'random') {
			({ vertices: newVertices, edges: newEdges } = generateRandomGraph(width, height));
		} else if (graphType === 'bipartite') {
			({ vertices: newVertices, edges: newEdges } = generate2x3BipartiteGraph(width, height));
		}
		$vertices = newVertices;
		$edges = newEdges;

		const initialKernel = initializeKernelState($vertices, $edges, alpha, beta);
		$kernelData = initialKernel;
		$previousEnergy = initialKernel.discreteEnergy;

		initialEdgeLengths = initialKernel.edgeProps.edgeLengths;

		optimizer = createOptimizer(
			$vertices,
			$edges,
			alpha,
			beta,
			initialKernel.disjointPairs,
			maxIterations,
			updateVisualization,
			initialEdgeLengths
		);

		cleanupInteractions();
		cleanupInteractions = setupInteractions(graphCanvas, $vertices, updateVisualization);

		// Reset transformations in the store
		canvasTransform.set({
			offsetX: 0,
			offsetY: 0,
			zoom: 1.0
		});
		updateVisualization();
	}

	function startOptimization() {
		if (optimizer && !isOptimizing) {
			console.log('Start optimization clicked');
			optimizer.start();
			isOptimizing = true;
		}
	}

	function stopOptimization() {
		if (optimizer && isOptimizing) {
			console.log('Stop optimization clicked');
			optimizer.stop();
			isOptimizing = false;
		}
	}

	function singleStep() {
		if (optimizer) {
			console.log('Single step clicked');
			optimizer.step();
		}
	}

	function updateAlphaBeta() {
		const updatedKernel = updateKernelState(
			$vertices,
			$edges,
			alpha,
			beta,
			$kernelData.disjointPairs
		);
		$kernelData = { ...updatedKernel, disjointPairs: $kernelData.disjointPairs };
		updateVisualization();
	}

	function updateConfig() {
		updateVisualization();
	}

	function getEnergyChangeColor() {
		return $energyChange < 0 ? 'green' : 'red';
	}
</script>

<div class="visualization-container">
	<div
		class="controls"
		style="display: flex; flex-direction: column; gap: 10px; top: 10px; left: 10px; z-index: 10;"
	>
		<button on:click={regenerateGraph}>Regenerate Graph</button>
		<button on:click={isOptimizing ? stopOptimization : startOptimization}>
			{isOptimizing ? 'Stop Optimization' : 'Start Optimization'}
		</button>
		<button on:click={singleStep}>Single Step</button>
		<div class="energy-value">
			<p>Discrete Energy: {$discreteEnergy.toFixed(4)}</p>
			<p style="color: {getEnergyChangeColor()}">Energy Change: {$energyChange.toFixed(4)}</p>
		</div>
		<Controls on:update={updateVisualization} />
	</div>
	<div class="graph-section">
		<div class="graph-container" style="position: relative; width: {width}px; height: {height}px;">
			<canvas bind:this={graphCanvas} {width} {height} style="position: absolute; top: 0; left: 0;"
			></canvas>
		</div>
	</div>
	<!-- <div class="kernel-section">
		<canvas bind:this={kernelCanvas}></canvas>
	</div> -->
</div>

<style>
	.visualization-container {
		display: flex;
		flex-direction: row;
		gap: 20px;
	}
	.graph-section,
	.kernel-section {
		flex: 1;
	}
</style>



---
File: /testdiscrete.js
---

import * as math from 'mathjs';

// Function to compute k_{β}^{α}(yi, yj, TI)
function k_alpha_beta(y_i, y_j, T_I, alpha, beta) {
	const diff = math.subtract(y_j, y_i);
	const crossProduct = math.norm(math.cross(T_I, diff));
	const distance = math.norm(diff);
	return Math.pow(crossProduct, alpha) / Math.pow(distance, beta);
}

// Compute the discrete tangent point energy
function computeDiscreteEnergy(edges, points, tangents, segmentLengths, alpha, beta) {
	let energy = 0;

	for (let I of edges) {
		for (let J of edges) {
			if (I !== J && I.filter((i) => J.includes(i)).length === 0) {
				// Ensure disjoint segments
				let sum_k = 0;

				for (let i of I) {
					for (let j of J) {
						sum_k += k_alpha_beta(points[i], points[j], tangents[I], alpha, beta);
					}
				}

				energy += (1 / 4) * sum_k * segmentLengths[I] * segmentLengths[J];
			}
		}
	}

	return energy;
}

// Example usage
const edges = [
	[0, 1],
	[2, 3]
]; // Example edge indices
const points = [
	[0, 0, 0],
	[1, 0, 0],
	[2, 1, 0],
	[3, 1, 0]
]; // Example 3D points
tangents = [
	[0, 1, 0],
	[0, 1, 0]
]; // Example tangents
const segmentLengths = { 0: 1, 1: 1 }; // Example segment lengths
const alpha = 2;
const beta = 3;

const energy = computeDiscreteEnergy(edges, points, tangents, segmentLengths, alpha, beta);
console.log('Discrete Tangent Point Energy:', energy);

---
---
let me remind you this important part of the paper:

Notation. In the discrete setting, we will model any collection of curves and loops (including several curves meeting at a common point) as a graph $G = (V,E)$ with vertex coordinates $\gamma : V \rightarrow \mathbb{R}^3$ (Figure 8); we use $|V|$ and $|E|$ to denote the number of vertices and edges, resp. For each edge $I \in E$ with endpoints $i_1, i_2$, we use

$\ell_I := |\gamma_{i_1} - \gamma_{i_2}|$, $T_I := (\gamma_{i_2} - \gamma_{i_1})/\ell_I$, and $x_I := (\gamma_{i_1} + \gamma_{i_2})/2$

to denote the edge length, unit tangent, and midpoint, resp. For any quantity $u: V \rightarrow \mathbb{R}$ on vertices we use $u_I := (u_{i_1} + u_{i_2})/2$ to denote the average value on edge $I = (i_1, i_2)$, and $u[I] := \begin{bmatrix} u_{i_1} & u_{i_2} \end{bmatrix}^T$ to denote the $2 \times 1$ column vector storing the values at its endpoints. Finally, we refer to any pair $(T,x) \in \mathbb{R}^6$ as a tangent-point.

5.1 Discrete Energy

Since the tangent-point energy is infinite for polygonal curves [Strzelecki and von der Mosel 2017, Figure 2.2], we assume that $\gamma$ is inscribed in some (unknown) smooth curve, and apply numerical quadrature to the smooth energy $\mathcal{E}_{\beta}^{\alpha}$. The resulting discrete energy then approximates the energy of any sufficiently smooth curve passing through the vertices $\gamma_i$. We start by integrating $k_{\beta}^{\alpha}$ over all pairs of edges:

$\sum_{I \in E} \sum_{J \in E} \int_{\overline{I}} \int_{\overline{J}} k_{\beta}^{\alpha}(\gamma(x), \gamma(y), T_I) dx_{\gamma} dy_{\gamma}.$ (16)

Here $\overline{I}$ denotes the interval along edge $I$. As given, this expression is ill-defined since two edges with a common endpoint contribute infinite energy. One idea is to instead use a term involving the curvature of the circle passing through the three distinct endpoints (in the spirit of Equation 1). However, such terms would contribute nothing to the energy in the limit of regular refinement (Figure 9) - hence, we simply omit neighboring edge pairs. Applying the (2D) trapezoidal rule to Equation 16 then yields a discrete energy

$\hat{\mathcal{E}}_{\beta}^{\alpha}(\gamma) = \sum_{I, J \in E, I \cap J = \emptyset} (\hat{k}_{\beta}^{\alpha})_{IJ} \ell_I \ell_J,$ (17)

where $\hat{k}$ is the discrete kernel

$(\hat{k}_{\beta}^{\alpha})_{IJ} := \frac{1}{4} \sum_{i \in I} \sum_{j \in J} k_{\beta}^{\alpha}(\gamma_i, \gamma_j, T_I)$. (18)

The discrete differential is then simply the partial derivatives of this energy with respect to the coordinates of all the curve vertices:

$$
d\hat{\mathcal{E}}^\alpha_\beta\big|_{\gamma} = \begin{bmatrix}
\partial \mathcal{E}^\alpha_\beta / \partial \gamma_1 & \cdots & \partial \mathcal{E}^\alpha_\beta / \partial \gamma_{|V|}
\end{bmatrix} \in \mathbb{R}^{3|V|}.
$$

These derivatives can be evaluated via any standard technique (*e.g.*, by hand, or using symbolic or automatic differentiation).

### 5.2 Discrete Inner Product

As in the smooth setting, we define our inner product matrix as a sum $A = B + B^0$ of high-order and low-order terms $B, B^0 \in \mathbb{R}^{|V| \times |V|}$ (as defined below).  For $\mathbb{R}^3$-valued functions, we also define a corresponding $3|V| \times 3|V|$ matrix

$$
\bar{A} =
\begin{bmatrix}
A & & \\
& A & \\
& & A
\end{bmatrix}. \qquad (19)
$$

Mirroring Equation 8, the discrete (fractional) Sobolev gradient $g \in \mathbb{R}^{3|V|}$ is then defined as the solution to the matrix equation

$$
\bar{A} g = d \hat{\mathcal{E}}_\beta^\alpha. \qquad (20)
$$

### 5.2.1 Discrete Derivative Operator.

For each edge \( I \in E \) we approximate the derivative \( \mathcal{D}u \) of a function \( u: M \to \mathbb{R} \) (Equation 11) via the finite difference formula \( \frac{1}{l_I}(u_{i_2} - u_{i_1})T_I \), where \( u_i \) denotes the value of \( u \) sampled at vertex \( i \). The corresponding derivative matrix \( D \in \mathbb{R}^{3|E| \times |V|} \) can be assembled from local \( 3 \times 2 \) matrices

\[
D_I = \frac{1}{l_I}
\begin{bmatrix}
-T_I & T_I
\end{bmatrix}.
\]

### 5.2.2 Discrete High-Order Term.

We approximate the high-order part of the inner product \( \langle \! \langle B_\sigma u, v \rangle \! \rangle \) as

\[
u^T B v = \sum_{I, J \in E, I \cap J = \emptyset} w_{IJ} \langle D_I u[I] - D_J u[J], D_I v[I] - D_J v[J] \rangle, \quad (21)
\]

where the weights \( w_{IJ} \) arise from applying trapezoidal quadrature to the denominator in Equation 25:

\[
w_{IJ} := \frac{1}{4} l_I l_J \sum_{i \in I} \sum_{j \in J} \frac{1}{|\gamma_i - \gamma_j|^{2\sigma + 1}}.
\]

The entries of the corresponding Gram matrix \( B \in \mathbb{R}^{|V| \times |V|} \) are obtained by differentiating Equation 21 with respect to the entries of \( u \) and \( v \). More explicitly, starting with the zero matrix one can build \( B \) by making the following increments for all pairs of disjoint edges \( I \cap J = \emptyset \), and all pairs of values \( a, b \in \{1, 2\} \):

\[
\begin{aligned}
B_{i_a i_b} += & (-1)^{a+b} w_{IJ} / l_I^2, & B_{i_a j_b} -= & (-1)^{a+b} w_{IJ} \langle T_I, T_J \rangle / (l_I l_J), \\
B_{j_a j_b} += & (-1)^{a+b} w_{IJ} / l_J^2, & B_{j_a i_b} -= & (-1)^{a+b} w_{IJ} \langle T_J, T_I \rangle / (l_J l_I). 
\end{aligned}
\]

### 5.2.3 Discrete Low-Order Term.

To discretize the low-order term \( B_\sigma^0 \) (Section 4.2.3), we use a different discrete weight

\[
w_{IJ}^0 := \frac{1}{4} l_I l_J \sum_{i \in I} \sum_{j \in J} \frac{k_4^2(\gamma_i, \gamma_j, T_I)}{|\gamma_i - \gamma_j|^{2\sigma + 1}},
\]

and define a matrix \( B^0 \in \mathbb{R}^{|V| \times |V|} \), given by the relationship

\[
u^T B^0 v = \sum_{I, J \in E, I \cap J = \emptyset} w_{IJ}^0 (u_I - u_J)(v_I - v_J).
\]

Following a similar derivation as above, this matrix can be constructed via the following increments:

\[
\begin{aligned}
B_{i_a i_b}^0 += & \frac{1}{4} w_{IJ}^0, & B_{i_a j_b}^0 -= & \frac{1}{4} w_{IJ}^0, \\
B_{j_a i_b}^0 -= & \frac{1}{4} w_{IJ}^0, & B_{j_a j_b}^0 += & \frac{1}{4} w_{IJ}^0.
\end{aligned}
\]
---
---
1. generate subvertices that are glued to edges and distribute them evenly across the edges
2. vary the amount of subvertices dynamically based on the length of the edge
2.1 vary the amount each step of the optimization and add a configuration how dense the vertices of edges should be i.e. their gap distance (by the config in the store). 
3. simple but Crucial: allow those vertices to be only inside the edge during optimization step without popping out, e.g. the vertex is allowed to be only the relative position inside the edge, cant move around in it, so that means only the supervertex(es) can move. i.e. if the subv. that is in the middle between two vertices is moved(the gradient is high) then kinda both superv. move. but if the subv is near some ver. then ofc that ver. is moved much more than the farther one. Perhaps it's a good idea to separate vertices with subverties to prevent them from moving freely because their behaviour is different. Debatable. This is all to make the energy minimization better while allowing to edges to interact with each other less. p.s. skip the interaction, i mostly care about the visual display of those vertices, forbid the free movement, and them contributing to the optimization. Provide the modified files in their entirety for easier copying except energycalculations.