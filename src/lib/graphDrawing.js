// src/lib/graphDrawing.js
import * as math from 'mathjs';
import { drawArrow } from '$lib/graphUtils';

export function drawGraph(ctx, width, height, vertices, edges, edgeProps, kernelMatrix) {
	ctx.clearRect(0, 0, width, height);

	// Draw edges with varying thickness and color based on kernel values
	if (kernelMatrix) {
		const maxKernelValue = math.max(kernelMatrix);
		console.log('Max kernel value:', maxKernelValue);

		for (let i = 0; i < edges.length; i++) {
			const edge = edges[i];

			// Calculate average kernel value for this edge
			let totalKernel = 0;
			let count = 0;
			for (let j = 0; j < edges.length; j++) {
				if (i !== j) {
					totalKernel += kernelMatrix.get([i, j]);
					count++;
				}
			}
			const avgKernel = totalKernel / (count || 1);
			const normalizedValue = avgKernel / maxKernelValue;

			console.log(`Edge ${i} avg kernel:`, avgKernel, 'normalized:', normalizedValue);

			// Calculate color (blue to red gradient)
			const blue = Math.round(255 * (1 - normalizedValue));
			const red = Math.round(255 * normalizedValue);

			// Set line properties
			ctx.strokeStyle = `rgb(${red}, 0, ${blue})`;
			ctx.lineWidth = 1 + normalizedValue * 4; // Thickness varies from 1 to 5 pixels

			// Draw edge
			ctx.beginPath();
			ctx.moveTo(vertices[edge[0]][0], vertices[edge[0]][1]);
			ctx.lineTo(vertices[edge[1]][0], vertices[edge[1]][1]);
			ctx.stroke();
		}
	}

	// Draw vertices
	for (let i = 0; i < vertices.length; i++) {
		// Draw vertex circle
		ctx.beginPath();
		ctx.arc(vertices[i][0], vertices[i][1], 5, 0, 2 * Math.PI);
		ctx.fillStyle = 'blue';
		ctx.fill();

		// Draw vertex number
		ctx.fillStyle = 'black';
		ctx.font = '12px Arial';
		ctx.textAlign = 'center';
		ctx.textBaseline = 'middle';
		ctx.fillText(i.toString(), vertices[i][0], vertices[i][1]);
	}

	// Draw midpoints with edge information and tangent arrows
	ctx.font = '10px Arial';
	ctx.textAlign = 'center';
	ctx.textBaseline = 'bottom';

	for (let i = 0; i < edges.length; i++) {
		const midpoint = edgeProps.edgeMidpoints[i];
		const length = edgeProps.edgeLengths[i].toFixed(2);
		const tangent = edgeProps.edgeTangents[i];

		// Draw midpoint dot
		ctx.fillStyle = 'red';
		ctx.beginPath();
		ctx.arc(midpoint[0], midpoint[1], 2, 0, 2 * Math.PI);
		ctx.fill();

		// Draw edge information
		ctx.fillStyle = 'blue';
		ctx.fillText(`L ${parseFloat(length).toFixed(0)}`, midpoint[0], midpoint[1] - 15);
		ctx.fillText(`T ${tangent.map((t) => t.toFixed(1))}`, midpoint[0], midpoint[1] - 5);
		ctx.fillText(`${edges[i][0]}, ${edges[i][1]}`, midpoint[0], midpoint[1] + 15);

		// Draw tangent arrow
		ctx.strokeStyle = 'green';
		ctx.lineWidth = 1.5;
		drawArrow(ctx, midpoint[0], midpoint[1], tangent[0], tangent[1], 20);
	}
}

export function drawKernelMatrix(kernelCtx, matrix) {
	if (!kernelCtx) return;

	const size = matrix.size();
	const padding = 50; // Padding for labels
	const maxWidth = Math.min(window.innerWidth / 2, 800); // Max width is half screen width or 800px
	const maxContentWidth = maxWidth - padding * 2; // Available space for cells

	// Calculate cell size to fit within maxWidth
	const cellSize = Math.min(
		50, // Maximum cell size
		Math.floor(maxContentWidth / size[1]), // Width-constrained cell size
		Math.floor(maxContentWidth / size[0]) // Height-constrained cell size (maintain square)
	);

	// Calculate actual dimensions
	const canvasWidth = Math.min(maxWidth, padding * 2 + size[1] * cellSize);
	const canvasHeight = padding * 2 + size[0] * cellSize;

	console.log('Matrix size:', size[0], 'x', size[1]);
	console.log('Cell size:', cellSize);
	console.log('Canvas dimensions:', canvasWidth, 'x', canvasHeight);

	// Update canvas dimensions
	kernelCtx.width = canvasWidth;
	kernelCtx.height = canvasHeight;

	kernelCtx.clearRect(0, 0, kernelCtx.width, kernelCtx.height);
	const maxVal = math.max(matrix);

	// Draw the matrix cells
	for (let i = 0; i < size[0]; i++) {
		for (let j = 0; j < size[1]; j++) {
			const value = matrix.get([i, j]);
			const intensity = value / maxVal;
			// Use white to blue gradient for better visibility
			const r = Math.round(255 - intensity * 255);
			const g = Math.round(255 - intensity * 255);
			const b = 255;
			const minOpacity = 0.1; // Ensure some minimum visibility
			const opacity = minOpacity + (1 - minOpacity) * intensity;
			const color = `rgba(${r}, ${g}, ${b}, ${opacity})`;

			kernelCtx.fillStyle = color;
			kernelCtx.fillRect(padding + j * cellSize, padding + i * cellSize, cellSize, cellSize);
		}
	}

	// Draw grid lines
	kernelCtx.strokeStyle = 'black';
	kernelCtx.lineWidth = 1;
	for (let i = 0; i <= size[0]; i++) {
		kernelCtx.beginPath();
		kernelCtx.moveTo(padding, padding + i * cellSize);
		kernelCtx.lineTo(padding + size[1] * cellSize, padding + i * cellSize);
		kernelCtx.stroke();
	}
	for (let j = 0; j <= size[1]; j++) {
		kernelCtx.beginPath();
		kernelCtx.moveTo(padding + j * cellSize, padding);
		kernelCtx.lineTo(padding + j * cellSize, padding + size[0] * cellSize);
		kernelCtx.stroke();
	}

	// Draw axis labels
	kernelCtx.fillStyle = 'black';
	kernelCtx.font = '12px Arial';
	kernelCtx.textAlign = 'center';

	// X axis labels (columns)
	for (let j = 0; j < size[1]; j++) {
		kernelCtx.fillText(j.toString(), padding + j * cellSize + cellSize / 2, padding - 10);
	}

	// Y axis labels (rows)
	kernelCtx.textAlign = 'right';
	for (let i = 0; i < size[0]; i++) {
		kernelCtx.fillText(i.toString(), padding - 5, padding + i * cellSize + cellSize / 2);
	}
}
