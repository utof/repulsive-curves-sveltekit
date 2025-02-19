<script>
	import { onMount, onDestroy } from 'svelte';
	import * as math from 'mathjs';
	import { calculateEdgeProperties, calculateDiscreteKernel } from '$lib/energyCalculations';

	let canvas;
	let ctx;
	let vertices = [];
	let edges = [];
	let edgeProps = { edgeLengths: [], edgeTangents: [], edgeMidpoints: [] };
	let kernelCanvas;
	let kernelCtx;

	const width = 700;
	const height = 700;
	const alpha = 3;
	const beta = 6;
	const TANGENT_ARROW_LENGTH = 20; // Length of tangent arrow visualization

	// Function to draw the kernel matrix visualization
	function drawKernelMatrix(matrix) {
		if (!kernelCtx) return;

		kernelCtx.clearRect(0, 0, kernelCanvas.width, kernelCanvas.height);
		const size = matrix.size();
		const maxVal = math.max(matrix);
		const cellSize = Math.min(kernelCanvas.width / size[0], kernelCanvas.height / size[1]);
		const padding = 30; // Space for axis labels

		// Draw the matrix cells
		for (let i = 0; i < size[0]; i++) {
			for (let j = 0; j < size[1]; j++) {
				const value = matrix.get([i, j]);
				const intensity = value / maxVal;
				const color = `rgba(0, 0, 255, ${intensity})`;

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

	onMount(() => {
		canvas = document.getElementById('graphCanvas');
		ctx = canvas.getContext('2d');
		kernelCanvas = document.getElementById('kernelCanvas');
		kernelCtx = kernelCanvas.getContext('2d');
		generateRandomGraph();
		calculateAndDrawKernel();
		drawGraph();
	});

	function generateRandomGraph() {
		const numVertices = Math.floor(Math.random() * 10) + 5;
		vertices = [];
		edges = [];

		for (let i = 0; i < numVertices; i++) {
			vertices.push([
				Math.random() * width, // x within canvas width
				Math.random() * height // y within canvas height
			]);
		}

		const edgeSet = new Set();
		for (let i = 0; i < numVertices; i++) {
			for (let j = i + 1; j < numVertices; j++) {
				if (Math.random() < 0.2) {
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
	}

	function calculateAndDrawKernel() {
		edgeProps = calculateEdgeProperties(vertices, edges);

		// Calculate and log discrete kernel
		const kernelMatrix = calculateDiscreteKernel(
			vertices,
			edges,
			edgeProps.edgeTangents,
			alpha,
			beta
		);
		console.log('Discrete Kernel Matrix:', kernelMatrix);
		drawKernelMatrix(kernelMatrix);
	}

	function drawArrow(fromX, fromY, dirX, dirY, length) {
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

	function drawGraph() {
		ctx.clearRect(0, 0, width, height);

		// Draw edges
		ctx.strokeStyle = 'gray';
		ctx.lineWidth = 1;
		for (const edge of edges) {
			ctx.beginPath();
			ctx.moveTo(vertices[edge[0]][0], vertices[edge[0]][1]);
			ctx.lineTo(vertices[edge[1]][0], vertices[edge[1]][1]);
			ctx.stroke();
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
			drawArrow(midpoint[0], midpoint[1], tangent[0], tangent[1], TANGENT_ARROW_LENGTH);
		}
	}

	// Animation functions commented out for now
	// function animate() {
	//   drawGraph();
	//   calculateAndDrawKernel();
	//   animationFrameId = requestAnimationFrame(animate);
	// }

	// function startAnimation() {
	//   animate();
	// }

	// function stopAnimation() {
	//   if (animationFrameId) {
	//     cancelAnimationFrame(animationFrameId);
	//   }
	// }

	function regenerateGraph() {
		generateRandomGraph();
		calculateAndDrawKernel();
		drawGraph();
	}
</script>

<div style="position: relative; width: {width}px; height: {height}px;">
	<canvas id="graphCanvas" {width} {height} style="position: absolute; top: 0; left: 0;"></canvas>
	<button on:click={regenerateGraph} style="position: absolute; top: 10px; left: 10px; z-index: 10;"
		>Regenerate Graph</button
	>
</div>

<canvas id="kernelCanvas" width="500" height="500"></canvas>
