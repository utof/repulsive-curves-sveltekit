<script>
	import Controls from '$lib/Controls.svelte';
	import { onMount, onDestroy } from 'svelte';
	import { drawGraph, drawKernelMatrix } from '$lib/graphDrawing';
	import { createOptimizer } from '$lib/optimization';
	import {
		generateRandomGraph,
		generateBipartiteGraph,
		generate2x3BipartiteGraph,
		generateCompleteSquareGraph,
		generateSimpleSquareGraph
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
		initialTotalLength,
		optimizationConfig
	} from '$lib/stores';
	import { calculateTotalLength } from '$lib/constraints';
	import { get } from 'svelte/store';

	let graphCanvas;
	let kernelCanvas;
	let graphCtx;
	let optimizer;
	let cleanupInteractions = () => {};
	let isOptimizing = false;
	let graphType = 'square';
	const width = 700;
	const height = 700;
	let controlsComponent;

	// Set initial date/time and username
	currentDateTimeUTC.set('2025-03-02 15:30:03');
	currentUser.set('utof');

	onMount(() => {
		graphCtx = graphCanvas.getContext('2d');
		regenerateGraph();
	});

	onDestroy(() => {
		if (optimizer) optimizer.stop();
		cleanupInteractions();
	});

	function updateVisualization() {
		// Get alpha and beta from the store
		const { alpha, beta } = get(optimizationConfig);

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

		if (kernelCanvas) {
			drawKernelMatrix(kernelCanvas, updatedKernel.kernelMatrix);
		}
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
		} else if (graphType === 'square') {
			({
				vertices: newVertices,
				edges: newEdges,
				subvertices: newSubvertices
			} = generateCompleteSquareGraph(width, height));
		} else if (graphType === 'simple') {
			({
				vertices: newVertices,
				edges: newEdges,
				subvertices: newSubvertices
			} = generateSimpleSquareGraph(width, height));
		}
		$vertices = newVertices;
		$edges = newEdges;
		$subvertices = newSubvertices;

		// Get alpha and beta from optimization config store
		const { alpha, beta } = get(optimizationConfig);

		const initialKernel = initializeKernelState($vertices, $edges, alpha, beta);
		$kernelData = initialKernel;
		$previousEnergy = initialKernel.discreteEnergy;

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
			1000, // maxIterations
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

	function handleAlphaBetaChange(event) {
		if (optimizer) {
			const { alpha, beta } = event.detail;
			optimizer.updateAlphaBeta(alpha, beta);
			updateVisualization();
		}
	}

	function getEnergyChangeColor() {
		return $energyChange < 0 ? 'green' : 'red';
	}

	function switchGraphType(type) {
		graphType = type;
		regenerateGraph();
	}
</script>

<div class="app-container">
	<header class="app-header">
		<h1>Tangent-Point Energy Optimization</h1>
		<div class="user-info"></div>
	</header>

	<div class="visualization-container">
		<div class="controls-section">
			<div class="action-buttons">
				<div class="graph-type-selector">
					<button
						class:active={graphType === 'bipartite'}
						on:click={() => switchGraphType('bipartite')}
					>
						Bipartite Graph
					</button>
					<button class:active={graphType === 'random'} on:click={() => switchGraphType('random')}>
						Random Graph
					</button>
				</div>

				<button on:click={regenerateGraph} class="main-button regenerate">
					<i class="icon-refresh"></i> Regenerate Graph
				</button>
				<button
					on:click={isOptimizing ? stopOptimization : startOptimization}
					class:active={isOptimizing}
					class="main-button {isOptimizing ? 'stop' : 'start'}"
				>
					<i class="icon-{isOptimizing ? 'stop' : 'play'}"></i>
					{isOptimizing ? 'Stop Optimization' : 'Start Optimization'}
				</button>
				<button on:click={singleStep} class="main-button step">
					<i class="icon-step"></i> Single Step
				</button>
			</div>

			<div class="energy-display">
				<h3>Energy Information</h3>
				<div class="energy-value">
					<label>Discrete Energy:</label>
					<span class="value">{$discreteEnergy.toFixed(4)}</span>
				</div>
				<div class="energy-value">
					<label>Energy Change:</label>
					<span class="value" style="color: {getEnergyChangeColor()}">
						{$energyChange.toFixed(4)}
					</span>
				</div>
			</div>

			<Controls
				on:update={updateVisualization}
				on:alphaBetaChange={handleAlphaBetaChange}
				bind:this={controlsComponent}
			/>
		</div>

		<div class="graph-section">
			<div class="graph-container">
				<canvas bind:this={graphCanvas} {width} {height} class="graph-canvas"></canvas>

				<div class="graph-controls">
					<div class="zoom-controls">
						<button title="Zoom In">+</button>
						<button title="Reset Zoom">⟳</button>
						<button title="Zoom Out">−</button>
					</div>
				</div>
			</div>

			<!-- Uncomment if you want to display kernel matrix visualization
			<div class="kernel-container">
				<h3>Kernel Matrix</h3>
				<canvas bind:this={kernelCanvas} width="200" height="200" class="kernel-canvas"></canvas>
			</div>
			-->
		</div>
	</div>
</div>

<style>
	.app-container {
		font-family:
			-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
		color: #333;
		max-width: 1200px;
		margin: 0 auto;
		padding: 20px;
	}

	.app-header {
		display: flex;
		justify-content: space-between;
		align-items: center;
		padding-bottom: 10px;
		margin-bottom: 20px;
		border-bottom: 1px solid #dee2e6;
	}

	.app-header h1 {
		margin: 0;
		font-size: 1.8rem;
		color: #2c3e50;
	}

	.user-info {
		display: flex;
		flex-direction: column;
		align-items: flex-end;
		font-size: 0.9rem;
		color: #6c757d;
	}

	.visualization-container {
		display: flex;
		gap: 20px;
	}

	.controls-section {
		width: 350px;
		flex-shrink: 0;
	}

	.action-buttons {
		margin-bottom: 15px;
	}

	.graph-type-selector {
		display: flex;
		margin-bottom: 10px;
	}

	.graph-type-selector button {
		flex: 1;
		padding: 8px;
		border: 1px solid #ced4da;
		background-color: #f8f9fa;
		cursor: pointer;
		font-size: 0.9rem;
	}

	.graph-type-selector button:first-child {
		border-top-left-radius: 4px;
		border-bottom-left-radius: 4px;
	}

	.graph-type-selector button:last-child {
		border-top-right-radius: 4px;
		border-bottom-right-radius: 4px;
	}

	.graph-type-selector button.active {
		background-color: #4caf50;
		color: white;
		border-color: #4caf50;
	}

	.main-button {
		display: block;
		width: 100%;
		padding: 10px;
		margin-bottom: 8px;
		border: none;
		border-radius: 4px;
		background-color: #6c757d;
		color: white;
		cursor: pointer;
		font-size: 1rem;
		transition: background-color 0.2s;
	}

	.main-button.regenerate {
		background-color: #6c757d;
	}

	.main-button.regenerate:hover {
		background-color: #5a6268;
	}

	.main-button.start {
		background-color: #28a745;
	}

	.main-button.start:hover {
		background-color: #218838;
	}

	.main-button.stop {
		background-color: #dc3545;
	}

	.main-button.stop:hover {
		background-color: #c82333;
	}

	.main-button.step {
		background-color: #007bff;
	}

	.main-button.step:hover {
		background-color: #0069d9;
	}

	.energy-display {
		background-color: #f8f9fa;
		border-radius: 8px;
		padding: 12px;
		margin-bottom: 15px;
		box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
	}

	.energy-display h3 {
		margin-top: 0;
		margin-bottom: 10px;
		font-size: 1.1rem;
		color: #495057;
	}

	.energy-value {
		display: flex;
		justify-content: space-between;
		margin-bottom: 5px;
		font-size: 1rem;
	}

	.energy-value label {
		font-weight: bold;
		color: #495057;
	}

	.energy-value .value {
		font-family: 'Courier New', monospace;
		font-weight: bold;
	}

	.graph-section {
		flex: 1;
		display: flex;
		flex-direction: column;
	}

	.graph-container {
		position: relative;
		width: 100%;
		height: 700px;
		background-color: #fff;
		border: 1px solid #dee2e6;
		border-radius: 8px;
		overflow: hidden;
		box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
	}

	.graph-canvas {
		position: absolute;
		top: 0;
		left: 0;
	}

	.graph-controls {
		position: absolute;
		top: 10px;
		right: 10px;
		z-index: 10;
	}

	.zoom-controls {
		display: flex;
		flex-direction: column;
		background: rgba(255, 255, 255, 0.8);
		border-radius: 4px;
		overflow: hidden;
		box-shadow: 0 1px 3px rgba(0, 0, 0, 0.2);
	}

	.zoom-controls button {
		border: none;
		background: none;
		font-size: 1.2rem;
		padding: 5px 10px;
		cursor: pointer;
		transition: background-color 0.2s;
	}

	.zoom-controls button:hover {
		background-color: #e9ecef;
	}

	.kernel-container {
		margin-top: 15px;
		padding: 15px;
		background-color: #fff;
		border: 1px solid #dee2e6;
		border-radius: 8px;
		box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
	}

	.kernel-container h3 {
		margin-top: 0;
		margin-bottom: 10px;
		font-size: 1.1rem;
		color: #495057;
	}

	.kernel-canvas {
		display: block;
		margin: 0 auto;
		max-width: 100%;
	}

	/* Basic icon styling - you would replace these with actual icons */
	.icon-refresh,
	.icon-play,
	.icon-stop,
	.icon-step {
		display: inline-block;
		width: 16px;
		height: 16px;
		margin-right: 5px;
		vertical-align: text-bottom;
	}

	@media (max-width: 900px) {
		.visualization-container {
			flex-direction: column;
		}

		.controls-section {
			width: 100%;
		}

		.graph-container {
			height: 500px;
		}
	}
</style>
