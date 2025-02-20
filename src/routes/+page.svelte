<script>
	import { onMount, onDestroy } from 'svelte';
	import * as math from 'mathjs';
	import {
		calculateEdgeProperties,
		calculateDisjointEdgePairs,
		calculateDiscreteKernel,
		calculateDiscreteEnergy
	} from '$lib/energyCalculations';

	let canvas;
	let ctx;
	let vertices = [];
	let edges = [];
	let edgeProps = { edgeLengths: [], edgeTangents: [], edgeMidpoints: [] };
	let kernelCanvas;
	let kernelCtx;
	let discreteEnergy = 0;
	let kernelMatrix = null;
	let disjointPairs = calculateDisjointEdgePairs(edges);

	const width = 700;
	const height = 700;
	// Higher alpha and lower beta for more pronounced differences
	const alpha = 3;
	const beta = 6;
	const TANGENT_ARROW_LENGTH = 20; // Length of tangent arrow visualization

	// Function to draw the kernel matrix visualization
	function drawKernelMatrix(matrix) {
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
		kernelCanvas.width = canvasWidth;
		kernelCanvas.height = canvasHeight;

		kernelCtx.clearRect(0, 0, kernelCanvas.width, kernelCanvas.height);
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

				// Add value label in cell for debugging
				kernelCtx.fillStyle = 'black';
				kernelCtx.font = '10px Arial';
				kernelCtx.textAlign = 'center';
				// kernelCtx.fillText(
				// 	value.toFixed(2),
				// 	padding + j * cellSize + cellSize / 2,
				// 	padding + i * cellSize + cellSize / 2
				// );
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
		disjointPairs = calculateDisjointEdgePairs(edges);
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
	}

	function calculateAndDrawKernel() {
		edgeProps = calculateEdgeProperties(vertices, edges);

		// Calculate and log discrete kernel and energy
		kernelMatrix = calculateDiscreteKernel(
			vertices,
			edges,
			edgeProps.edgeTangents,
			alpha,
			beta,
			disjointPairs
		);
		discreteEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);

		console.log('Discrete Kernel Matrix:', kernelMatrix);
		console.log('Discrete Energy:', discreteEnergy);
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
		disjointPairs = calculateDisjointEdgePairs(edges); // Move this after graph generation
		calculateAndDrawKernel();
		drawGraph();
	}
</script>

<div class="visualization-container">
	<div class="graph-section">
		<div class="graph-container" style="position: relative; width: {width}px; height: {height}px;">
			<div class="energy-value" style="position: absolute; top: 10px; left: 50px; z-index: 10;">
				<p>Discrete Energy: {discreteEnergy.toFixed(4)}</p>
			</div>
			<canvas id="graphCanvas" {width} {height} style="position: absolute; top: 0; left: 0;"
			></canvas>
			<button
				on:click={regenerateGraph}
				style="position: absolute; top: 10px; left: 10px; z-index: 10;">Regenerate Graph</button
			>
		</div>
	</div>

	<div class="kernel-section">
		<canvas id="kernelCanvas"></canvas>
	</div>
</div>

<style>
	.visualization-container {
		display: flex;
		gap: 20px;
		align-items: flex-start;
		padding: 20px;
		max-width: 100%;
		overflow-x: auto;
	}

	.graph-section {
		flex-shrink: 0;
	}

	.kernel-section {
		flex-shrink: 0;
		display: flex;
		justify-content: center;
	}

	canvas#kernelCanvas {
		display: block;
		background: white;
	}

	.energy-value {
		font-size: 16px;
		font-weight: bold;
		text-align: center;
	}
</style>
