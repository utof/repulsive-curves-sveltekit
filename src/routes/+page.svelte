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
	const stepSize = 1000;
	const maxIterations = 1000;
	let optimizer = null;
	let cleanupInteractions;

	onMount(() => {
		canvas = document.getElementById('graphCanvas');
		ctx = canvas.getContext('2d');
		kernelCanvas = document.getElementById('kernelCanvas');

		regenerateGraph(); // Initial graph

		cleanupInteractions = setupInteractions(canvas, vertices, updateVisualization, width, height);
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
			// Call kernel.update() and *destructure* the results
			const { kernelMatrix, discreteEnergy: newDiscreteEnergy } = kernel.update();
			energyChange = newDiscreteEnergy - previousEnergy;
			previousEnergy = newDiscreteEnergy;
			discreteEnergy = newDiscreteEnergy; // Update the displayed energy

			drawGraph(
				ctx,
				width,
				height,
				vertices,
				edges,
				kernel.edgeProps, // Access edgeProps directly
				kernelMatrix // Use the updated kernelMatrix
			);
		} else {
			// Handle the case where kernel is not yet initialized (shouldn't happen normally)
			drawGraph(
				ctx,
				width,
				height,
				vertices,
				edges,
				{ edgeLengths: [], edgeTangents: [], edgeMidpoints: [] },
				null
			);
		}
	}

	function regenerateGraph() {
		previousEnergy = 0;
		energyChange = 0;
		const newGraph = generateRandomGraph(width, height);
		vertices = newGraph.vertices;
		edges = newGraph.edges;

		// Initialize or update kernel.  Crucially, *always* call setupKernel.
		kernel = setupKernel(kernelCanvas, vertices, edges, alpha, beta);

		// Create/Recreate the optimizer.  This is important because the optimizer
		// needs to work with the *current* vertices and edges.
		if (optimizer) {
			optimizer.stop(); // Stop any existing optimizer
		}
		optimizer = createOptimizer(
			vertices,
			edges,
			alpha,
			beta,
			kernel.disjointPairs,
			stepSize,
			width,
			height,
			maxIterations,
			updateVisualization
		);

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
