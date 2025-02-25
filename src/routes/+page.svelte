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
