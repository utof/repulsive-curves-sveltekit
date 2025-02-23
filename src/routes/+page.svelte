<script>
	import { onMount, onDestroy } from 'svelte';
	import {
		calculateEdgeProperties,
		calculateDisjointEdgePairs,
		calculateDiscreteKernel,
		calculateDiscreteEnergy,
		calculateL2Gradient
	} from '$lib/energyCalculations';
	import { drawGraph, drawKernelMatrix } from '$lib/graphDrawing';
	import { createOptimizer } from '$lib/optimization';
	import { generateRandomGraph, drawArrow } from '$lib/graphUtils';

	let canvas;
	let ctx;
	let vertices = [];
	let edges = [];
	let edgeProps = { edgeLengths: [], edgeTangents: [], edgeMidpoints: [] };
	let kernelCanvas;
	let kernelCtx;
	let discreteEnergy = 0;
	let previousEnergy = 0; // Store the previous energy value
	let energyChange = 0;
	let kernelMatrix = null;
	let disjointPairs = [];

	const width = 700;
	const height = 700;
	const alpha = 3;
	const beta = 6;
	const stepSize = 0.01;
	const maxIterations = 100;
	let optimizer = null;

	// Dragging state
	let draggingVertex = null;
	let dragOffsetX = 0;
	let dragOffsetY = 0;

	onMount(() => {
		canvas = document.getElementById('graphCanvas');
		ctx = canvas.getContext('2d');
		kernelCanvas = document.getElementById('kernelCanvas');
		kernelCtx = kernelCanvas.getContext('2d');

		regenerateGraph(); // Initial graph

		// Add event listeners for mouse events
		canvas.addEventListener('mousedown', handleMouseDown);
		canvas.addEventListener('mousemove', handleMouseMove);
		canvas.addEventListener('mouseup', handleMouseUp);
		canvas.addEventListener('mouseleave', handleMouseUp); // Stop dragging if mouse leaves canvas

		optimizer = createOptimizer(
			vertices,
			edges,
			alpha,
			beta,
			disjointPairs,
			stepSize,
			width,
			height,
			maxIterations,
			() => {
				// onUpdate callback
				disjointPairs = calculateDisjointEdgePairs(edges);
				calculateAndDrawKernel();
				drawGraph(ctx, width, height, vertices, edges, edgeProps, kernelMatrix);
			}
		);
	});

	onDestroy(() => {
		if (optimizer) {
			optimizer.stop(); // Clean up the optimizer
		}
		// Remove event listeners
		if (canvas) {
			canvas.removeEventListener('mousedown', handleMouseDown);
			canvas.removeEventListener('mousemove', handleMouseMove);
			canvas.removeEventListener('mouseup', handleMouseUp);
			canvas.removeEventListener('mouseleave', handleMouseUp);
		}
	});

	function calculateAndDrawKernel() {
		previousEnergy = discreteEnergy; // Store before recalculating
		edgeProps = calculateEdgeProperties(vertices, edges);
		kernelMatrix = calculateDiscreteKernel(
			vertices,
			edges,
			edgeProps.edgeTangents,
			alpha,
			beta,
			disjointPairs
		);
		discreteEnergy = calculateDiscreteEnergy(vertices, edges, alpha, beta, disjointPairs);
		energyChange = discreteEnergy - previousEnergy; // Calculate the change
		drawKernelMatrix(kernelCtx, kernelMatrix);
	}

	function regenerateGraph() {
		previousEnergy = 0; // Reset previous energy
		energyChange = 0;
		const newGraph = generateRandomGraph(width, height);
		vertices = newGraph.vertices;
		edges = newGraph.edges;
		disjointPairs = calculateDisjointEdgePairs(edges);
		calculateAndDrawKernel();
		drawGraph(ctx, width, height, vertices, edges, edgeProps, kernelMatrix);
	}

	function startOptimization() {
		if (optimizer) {
			optimizer.start();
		}
	}

	function singleStep() {
		if (optimizer) {
			optimizer.step();
		}
	}

	// Mouse event handlers
	function handleMouseDown(event) {
		const rect = canvas.getBoundingClientRect();
		const mouseX = event.clientX - rect.left;
		const mouseY = event.clientY - rect.top;

		// Check if the mouse is over a vertex
		for (let i = 0; i < vertices.length; i++) {
			const [vx, vy] = vertices[i];
			const distance = Math.sqrt((mouseX - vx) ** 2 + (mouseY - vy) ** 2);
			if (distance <= 5) {
				// 5 is the vertex radius
				draggingVertex = i;
				dragOffsetX = vx - mouseX;
				dragOffsetY = vy - mouseY;
				break;
			}
		}
	}

	function handleMouseMove(event) {
		if (draggingVertex !== null) {
			const rect = canvas.getBoundingClientRect();
			const mouseX = event.clientX - rect.left;
			const mouseY = event.clientY - rect.top;

			// Update vertex position, constrained to canvas bounds
			let newX = mouseX + dragOffsetX;
			let newY = mouseY + dragOffsetY;

			newX = Math.max(0, Math.min(width, newX)); //Clamp inside
			newY = Math.max(0, Math.min(height, newY));

			vertices[draggingVertex] = [newX, newY];

			// Recalculate and redraw
			calculateAndDrawKernel();
			drawGraph(ctx, width, height, vertices, edges, edgeProps, kernelMatrix);
		}
	}

	function handleMouseUp(event) {
		draggingVertex = null;
	}

	function getEnergyChangeColor() {
		return energyChange < 0 ? 'green' : 'red';
	}
</script>

<div class="visualization-container">
	<div class="graph-section">
		<div class="graph-container" style="position: relative; width: {width}px; height: {height}px;">
			<div
				class="controls"
				style="display: flex; flex-direction: column; gap: 10px; position: absolute; top: 10px; left: 10px; z-index: 10;"
			>
				<button on:click={regenerateGraph}>Regenerate Graph</button>
				<button on:click={startOptimization}>Start Optimization</button>
				<button on:click={singleStep}>Single Step</button>
				<div class="energy-value">
					<p>Discrete Energy: {discreteEnergy.toFixed(4)}</p>
					<p style="color: {getEnergyChangeColor()}">
						Energy Change: {energyChange.toFixed(4)}
					</p>
				</div>
			</div>
			<canvas id="graphCanvas" {width} {height} style="position: absolute; top: 0; left: 0;"
			></canvas>
		</div>
	</div>

	<div class="kernel-section">
		<canvas id="kernelCanvas"></canvas>
	</div>
	content_copy download
</div>

<style>
	.visualization-container {
		display: flex;
		flex-direction: row; /* or column, depending on layout */
		gap: 20px;
	}

	.graph-section,
	.kernel-section {
		flex: 1; /* Allow both sections to grow and take available space */
	}

	/* If you need a fixed width for the kernel section, you can do this: */
	/* .kernel-section { */
	/*   width: 400px;  Adjust as needed */
	/* } */
</style>
