<script>
	import { onMount, onDestroy } from 'svelte';
	import { drawGraph, setupKernel } from '$lib/graphDrawing'; // Single-line import
	import { createOptimizer } from '$lib/optimization';
	import { generateRandomGraph } from '$lib/graphUtils';
	import { setupInteractions } from '$lib/interaction';

	let canvas;
	let ctx;
	let vertices = [];
	let edges = [];
	let kernelCanvas;
	let kernel; // Store kernel object
	let discreteEnergy = 0;
	let previousEnergy = 0;
	let energyChange = 0;

	const width = 700;
	const height = 700;
	const alpha = 3;
	const beta = 6;
	const stepSize = 0.01;
	const maxIterations = 100;
	let optimizer = null;
	let cleanupInteractions;

	onMount(() => {
		canvas = document.getElementById('graphCanvas');
		ctx = canvas.getContext('2d');
		kernelCanvas = document.getElementById('kernelCanvas');

		regenerateGraph(); // Initial graph

		cleanupInteractions = setupInteractions(canvas, vertices, updateVisualization, width, height);

		optimizer = createOptimizer(
			vertices,
			edges,
			alpha,
			beta,
			kernel ? kernel.disjointPairs : [], // Use kernel.disjointPairs if available
			stepSize,
			width,
			height,
			maxIterations,
			updateVisualization // Use a single update function
		);
	});

	onDestroy(() => {
		if (optimizer) {
			optimizer.stop();
		}
		if (cleanupInteractions) {
			cleanupInteractions();
		}
	});

	function updateVisualization() {
		if (kernel) {
			const updatedKernelData = kernel.update(); // Call the update function
			discreteEnergy = updatedKernelData.discreteEnergy;
			energyChange = discreteEnergy - previousEnergy;
			previousEnergy = discreteEnergy;
		}
		drawGraph(
			ctx,
			width,
			height,
			vertices,
			edges,
			kernel ? kernel.edgeProps : { edgeLengths: [], edgeTangents: [], edgeMidpoints: [] },
			kernel ? kernel.kernelMatrix : null
		);
	}

	function regenerateGraph() {
		previousEnergy = 0;
		energyChange = 0;
		const newGraph = generateRandomGraph(width, height);
		vertices = newGraph.vertices;
		edges = newGraph.edges;

		// Initialize or update kernel
		kernel = kernel ? kernel.update() : setupKernel(kernelCanvas, vertices, edges, alpha, beta);
		updateVisualization(); // Initial draw
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
</style>
