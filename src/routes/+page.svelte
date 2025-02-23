<!-- src/routes/+page.svelte -->
<script>
	import { onMount, onDestroy } from 'svelte';
	import { drawGraph, drawKernelMatrix } from '$lib/graphDrawing';
	import { createOptimizer } from '$lib/optimization';
	import { generateRandomGraph } from '$lib/graphUtils';
	import { setupInteractions } from '$lib/interaction';
	import { initializeKernelState, updateKernelState } from '$lib/graphState';
	import {
		vertices,
		edges,
		kernelData,
		energyChange,
		previousEnergy,
		discreteEnergy
	} from '$lib/stores';

	let graphCanvas;
	let kernelCanvas;
	let graphCtx;
	let optimizer;
	let cleanupInteractions;

	const width = 700;
	const height = 700;
	const alpha = 3;
	const beta = 6;
	const stepSize = 1000;
	const maxIterations = 1000;

	onMount(() => {
		graphCtx = graphCanvas.getContext('2d');
		regenerateGraph();
		cleanupInteractions = setupInteractions(
			graphCanvas,
			$vertices,
			updateVisualization,
			width,
			height
		);
	});

	onDestroy(() => {
		if (optimizer) optimizer.stop();
		if (cleanupInteractions) cleanupInteractions();
	});

	function updateVisualization() {
		const updatedKernel = updateKernelState(
			$vertices,
			$edges,
			alpha,
			beta,
			$kernelData.disjointPairs
		);
		kernelData.set({ ...updatedKernel, disjointPairs: $kernelData.disjointPairs });

		$energyChange = $discreteEnergy - $previousEnergy;
		$previousEnergy = $discreteEnergy;

		drawGraph(
			graphCtx,
			width,
			height,
			$vertices,
			$edges,
			updatedKernel.edgeProps,
			updatedKernel.kernelMatrix
		);
		drawKernelMatrix(kernelCanvas, updatedKernel.kernelMatrix);
	}

	function regenerateGraph() {
		$previousEnergy = 0;
		$energyChange = 0;
		const { vertices: newVertices, edges: newEdges } = generateRandomGraph(width, height);
		$vertices = newVertices;
		$edges = newEdges;

		const initialKernel = initializeKernelState($vertices, $edges, alpha, beta);
		$kernelData = initialKernel;
		$previousEnergy = initialKernel.discreteEnergy;

		if (optimizer) optimizer.stop();
		optimizer = createOptimizer(
			$vertices,
			$edges,
			alpha,
			beta,
			initialKernel.disjointPairs,
			stepSize,
			width,
			height,
			maxIterations,
			updateVisualization
		);

		updateVisualization();
	}

	function startOptimization() {
		if (optimizer) optimizer.start();
	}

	function singleStep() {
		if (optimizer) optimizer.step();
	}

	function getEnergyChangeColor() {
		return $energyChange < 0 ? 'green' : 'red';
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
					<p>Discrete Energy: {$discreteEnergy.toFixed(4)}</p>
					<p style="color: {getEnergyChangeColor()}">Energy Change: {$energyChange.toFixed(4)}</p>
				</div>
			</div>
			<canvas bind:this={graphCanvas} {width} {height} style="position: absolute; top: 0; left: 0;"
			></canvas>
		</div>
	</div>
	<div class="kernel-section">
		<canvas bind:this={kernelCanvas}></canvas>
	</div>
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
