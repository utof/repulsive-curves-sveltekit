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

function drawVertices(ctx, vertices, edges, alpha, beta, disjointPairs) {
	const { offsetX, offsetY, zoom } = get(canvasTransform);
	const gradient = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);

	vertices.forEach((vertex, i) => {
		// Draw circle in world coordinates
		ctx.beginPath();
		ctx.arc(vertex[0], vertex[1], 5 / zoom, 0, 2 * Math.PI); // World coordinates, scaled for zoom
		ctx.fillStyle = 'blue';
		ctx.fill();

		// Transform to screen coordinates for text
		const screenX = vertex[0] * zoom + offsetX;
		const screenY = vertex[1] * zoom + offsetY;
		ctx.fillStyle = 'black';
		ctx.font = '12px Arial';
		ctx.textAlign = 'center';
		ctx.textBaseline = 'bottom';
		ctx.fillText(i.toString(), screenX, screenY - 10 / zoom); // Adjust for zoom

		// Draw gradient arrow in world coordinates, transform to screen
		const gradX = -gradient[i][0]; // Negate for -gradient direction
		const gradY = -gradient[i][1];
		const magnitude = Math.sqrt(gradX * gradX + gradY * gradY) || 1e-6; // Avoid division by zero
		const minLength = 20; // Min arrow length for normalization
		const maxLength = 50; // Max arrow length for normalization
		const normalizedLength = Math.max(minLength, Math.min(maxLength, magnitude * 1000)) / zoom; // Scale for zoom

		const arrowStartX = vertex[0];
		const arrowStartY = vertex[1];
		const arrowEndX = arrowStartX + (gradX / magnitude) * normalizedLength;
		const arrowEndY = arrowStartY + (gradY / magnitude) * normalizedLength;

		// Transform arrow endpoints to screen coordinates
		const screenStartX = arrowStartX * zoom + offsetX;
		const screenStartY = arrowStartY * zoom + offsetY;
		const screenEndX = arrowEndX * zoom + offsetX;
		const screenEndY = arrowEndY * zoom + offsetY;

		ctx.strokeStyle = 'purple'; // Distinct color for gradient arrows
		ctx.lineWidth = 2 / zoom;
		drawArrow(
			ctx,
			screenStartX,
			screenStartY,
			(screenEndX - screenStartX) / normalizedLength,
			(screenEndY - screenStartY) / normalizedLength,
			normalizedLength
		);
	});
}

function drawMidpoints(ctx, edges, edgeProps) {
	const { offsetX, offsetY, zoom } = get(canvasTransform);

	if (!edgeProps || !edgeProps.edgeMidpoints || edgeProps.edgeMidpoints.length !== edges.length) {
		return;
	}

	ctx.font = '10px Arial';
	ctx.textAlign = 'center';
	ctx.textBaseline = 'bottom';

	edges.forEach((edge, i) => {
		const midpoint = edgeProps.edgeMidpoints[i];
		const length = edgeProps.edgeLengths[i];
		const tangent = edgeProps.edgeTangents[i];

		if (!midpoint || !length || !tangent) return;

		// Draw midpoint circle in world coordinates
		ctx.fillStyle = 'red';
		ctx.beginPath();
		ctx.arc(midpoint[0], midpoint[1], 2 / zoom, 0, 2 * Math.PI); // World coordinates, scaled for zoom
		ctx.fill();

		// Transform to screen coordinates for text and arrows
		const screenX = midpoint[0] * zoom + offsetX;
		const screenY = midpoint[1] * zoom + offsetY;

		ctx.fillStyle = 'blue';
		ctx.fillText(`${i}L ${length.toFixed(6)}`, screenX, screenY - 15 / zoom); // Adjust for zoom
		ctx.fillText(`T ${tangent.map((t) => t.toFixed(1))}`, screenX, screenY - 5 / zoom);
		ctx.fillText(`${edge[0]}, ${edge[1]}`, screenX, screenY + 15 / zoom);

		ctx.strokeStyle = 'green';
		ctx.lineWidth = 1.5 / zoom;
		drawArrow(ctx, screenX, screenY, tangent[0], tangent[1], 20 / zoom); // Adjust for zoom
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
