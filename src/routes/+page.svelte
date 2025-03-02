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
		subvertices,
		kernelData,
		energyChange,
		previousEnergy,
		discreteEnergy,
		config,
		canvasTransform,
		initialTotalLength
	} from '$lib/stores';
	import { calculateTotalLength } from '$lib/constraints';
	import { get } from 'svelte/store';

	let graphCanvas;
	let kernelCanvas;
	let graphCtx;
	let optimizer;
	let cleanupInteractions = () => {};
	let isOptimizing = false;
	let graphType = 'bipartite';
	const width = 700;
	const height = 700;
	let alpha = 3;
	let beta = 6;
	const maxIterations = 1000;
	let initialEdgeLengths = [];
	let controlsComponent;

	onMount(() => {
		graphCtx = graphCanvas.getContext('2d');
		regenerateGraph();
	});

	onDestroy(() => {
		if (optimizer) optimizer.stop();
		cleanupInteractions();
	});

	function handleAlphaBetaChange(event) {
		alpha = event.detail.alpha;
		beta = event.detail.beta;
		updateAlphaBeta();
	}

	function updateVisualization() {
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

		drawKernelMatrix(kernelCanvas, updatedKernel.kernelMatrix);
	}

	function regenerateGraph() {
		$previousEnergy = 0;
		$energyChange = 0;
		let newVertices, newEdges, newSubvertices;
		if (graphType === 'random') {
			({
				vertices: newVertices,
				edges: newEdges,
				subvertices: newSubvertices
			} = generateRandomGraph(width, height));
		} else if (graphType === 'bipartite') {
			({
				vertices: newVertices,
				edges: newEdges,
				subvertices: newSubvertices
			} = generate2x3BipartiteGraph(width, height));
		}
		$vertices = newVertices;
		$edges = newEdges;
		$subvertices = newSubvertices;

		const initialKernel = initializeKernelState($vertices, $edges, alpha, beta);
		$kernelData = initialKernel;
		$previousEnergy = initialKernel.discreteEnergy;

		initialEdgeLengths = initialKernel.edgeProps.edgeLengths;

		// Calculate and store the initial total curve length for percentage-based constraints
		const initTotalLength = calculateTotalLength($vertices, $edges);
		initialTotalLength.set(initTotalLength);
		console.log(`Initial total curve length: ${initTotalLength}`);

		optimizer = createOptimizer(
			$vertices,
			$edges,
			alpha,
			beta,
			initialKernel.disjointPairs,
			maxIterations,
			updateVisualization
		);

		// Connect the optimizer to the controls
		if (controlsComponent) {
			controlsComponent.setOptimizer(optimizer);
		}

		cleanupInteractions();
		cleanupInteractions = setupInteractions(graphCanvas, $vertices, updateVisualization);

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
		<Controls
			on:update={updateVisualization}
			on:alphaBetaChange={handleAlphaBetaChange}
			bind:this={controlsComponent}
		/>
	</div>
	<div class="graph-section">
		<div class="graph-container" style="position: relative; width: {width}px; height: {height}px;">
			<canvas bind:this={graphCanvas} {width} {height} style="position: absolute; top: 0; left: 0;"
			></canvas>
		</div>
	</div>
</div>

<style>
	.visualization-container {
		display: flex;
		flex-direction: row;
		gap: 20px;
	}
	.graph-section {
		flex: 1;
	}

	.energy-value {
		background-color: #f5f5f5;
		padding: 10px;
		border-radius: 5px;
		margin-bottom: 10px;
	}

	.energy-value p {
		margin: 0;
		padding: 3px 0;
	}

	button {
		padding: 8px 16px;
		font-size: 14px;
		border: none;
		border-radius: 4px;
		background-color: #4caf50;
		color: white;
		cursor: pointer;
		transition: background-color 0.3s;
	}

	button:hover {
		background-color: #45a049;
	}

	button:active {
		background-color: #3e8e41;
	}
</style>
