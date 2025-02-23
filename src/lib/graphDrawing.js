// src/lib/graphDrawing.js
import * as math from 'mathjs';
import { drawArrow } from '$lib/graphUtils';
import { calculateL2Gradient } from '$lib/energyCalculations';

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
	// Added parameters
	ctx.clearRect(0, 0, width, height);

	drawEdges(ctx, vertices, edges, kernelMatrix);
	drawVertices(ctx, vertices, edges, alpha, beta, disjointPairs); // Updated call
	drawMidpoints(ctx, edges, edgeProps);
}

function drawEdges(ctx, vertices, edges, kernelMatrix) {
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
			ctx.lineWidth = 1 + normalizedValue * 4;

			ctx.beginPath();
			ctx.moveTo(vertices[edge[0]][0], vertices[edge[0]][1]);
			ctx.lineTo(vertices[edge[1]][0], vertices[edge[1]][1]);
			ctx.stroke();
		});
	} else {
		// Fallback: Draw edges with default style
		ctx.strokeStyle = 'black';
		ctx.lineWidth = 1;
		edges.forEach((edge) => {
			ctx.beginPath();
			ctx.moveTo(vertices[edge[0]][0], vertices[edge[0]][1]);
			ctx.lineTo(vertices[edge[1]][0], vertices[edge[1]][1]);
			ctx.stroke();
		});
	}
}

function drawVertices(ctx, vertices, edges, alpha, beta, disjointPairs) {
	// Added parameters
	const gradient = calculateL2Gradient(vertices, edges, alpha, beta, disjointPairs);

	vertices.forEach((vertex, i) => {
		// Draw circle
		ctx.beginPath();
		ctx.arc(vertex[0], vertex[1], 5, 0, 2 * Math.PI);
		ctx.fillStyle = 'blue';
		ctx.fill();

		// Draw vertex index above
		ctx.fillStyle = 'black';
		ctx.font = '12px Arial';
		ctx.textAlign = 'center';
		ctx.textBaseline = 'bottom';
		ctx.fillText(i.toString(), vertex[0], vertex[1] - 10);

		// Draw gradient arrow
		const gradX = -gradient[i][0]; // Negate for -gradient direction
		const gradY = -gradient[i][1];
		const magnitude = Math.sqrt(gradX * gradX + gradY * gradY) || 1e-6; // Avoid division by zero
		const minLength = 20; // Min arrow length for normalization
		const maxLength = 50; // Max arrow length for normalization
		const normalizedLength = Math.max(minLength, Math.min(maxLength, magnitude * 1000)); // Scale factor (tweakable)

		ctx.strokeStyle = 'purple'; // Distinct color for gradient arrows
		ctx.lineWidth = 2;
		drawArrow(ctx, vertex[0], vertex[1], gradX / magnitude, gradY / magnitude, normalizedLength);
	});
}

function drawMidpoints(ctx, edges, edgeProps) {
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

		ctx.fillStyle = 'red';
		ctx.beginPath();
		ctx.arc(midpoint[0], midpoint[1], 2, 0, 2 * Math.PI);
		ctx.fill();

		ctx.fillStyle = 'blue';
		ctx.fillText(`${i}L ${length.toFixed(6)}`, midpoint[0], midpoint[1] - 15); // Modified line
		ctx.fillText(`T ${tangent.map((t) => t.toFixed(1))}`, midpoint[0], midpoint[1] - 5);
		ctx.fillText(`${edge[0]}, ${edge[1]}`, midpoint[0], midpoint[1] + 15);

		ctx.strokeStyle = 'green';
		ctx.lineWidth = 1.5;
		drawArrow(ctx, midpoint[0], midpoint[1], tangent[0], tangent[1], 20);
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
