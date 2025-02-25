<script>
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
		canvasTransform
	} from '$lib/stores';
	import { get } from 'svelte/store';
	import * as THREE from 'three';
	import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls';

	let graphCanvas2D;
	let graphCanvas3D;
	let graphCtx2D;
	let threeRenderer;
	let threeScene;
	let threeCamera;
	let threeControls;
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

	onMount(() => {
		vertices.set([]);
		edges.set([]);
		subvertices.set([]);
		kernelData.set({
			kernelMatrix: null,
			discreteEnergy: 0,
			edgeProps: { edgeLengths: [], edgeTangents: [], edgeMidpoints: [] },
			disjointPairs: []
		});
		energyChange.set(0);
		previousEnergy.set(0);
		config.set({
			dim: 2,
			epsilonStability: 1e-7,
			epsilonKernel: 1e-6,
			finiteDiffH: 1e-4,
			constraintTolerance: 1e-7,
			tauInitial: 1.0,
			aConst: 0.1,
			bConst: 0.5,
			maxLineSearch: 20,
			differentialMethod: 'finiteDifference',
			precondStepSize: 20,
			l2StepSize: 100000,
			applyPerturbation: false,
			subvertexGap: 50
		});
		canvasTransform.set({ offsetX: 0, offsetY: 0, zoom: 1.0 });

		graphCtx2D = graphCanvas2D.getContext('2d');
		setupThreeJS();
		regenerateGraph();

		function animate() {
			requestAnimationFrame(animate);
			if ($config.dim === 3) {
				updateVisualization();
			}
		}
		animate();
	});

	onDestroy(() => {
		if (optimizer) optimizer.stop();
		cleanupInteractions();
		if (threeRenderer) threeRenderer.dispose();
	});

	function setupThreeJS() {
		threeScene = new THREE.Scene();
		// Adjusted near and far clipping planes
		threeCamera = new THREE.PerspectiveCamera(75, width / height, 0.01, 10000); // Near: 0.01, Far: 10000
		threeCamera.position.set(0, 0, 500);
		threeRenderer = new THREE.WebGLRenderer({ canvas: graphCanvas3D });
		threeRenderer.setSize(width, height);
		threeControls = new OrbitControls(threeCamera, threeRenderer.domElement);
		threeControls.enableDamping = true;
		// Set min/max distance to prevent zooming too far in/out
		threeControls.minDistance = 50; // Minimum zoom distance
		threeControls.maxDistance = 2000; // Maximum zoom distance
		threeControls.update();
	}

	function updateVisualization() {
		const updatedKernel = updateKernelState(
			$vertices,
			$edges,
			alpha,
			beta,
			$kernelData.disjointPairs
		);
		kernelData.set({ ...updatedKernel, disjointPairs: $kernelData.disjointPairs });
		energyChange.set($discreteEnergy - $previousEnergy);
		previousEnergy.set($discreteEnergy);

		drawGraph(
			graphCtx2D,
			width,
			height,
			$vertices,
			$edges,
			updatedKernel.edgeProps,
			updatedKernel.kernelMatrix,
			alpha,
			beta,
			$kernelData.disjointPairs,
			threeRenderer,
			threeScene,
			threeCamera,
			threeControls
		);
	}

	function regenerateGraph() {
		previousEnergy.set(0);
		energyChange.set(0);
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
		vertices.set(newVertices);
		edges.set(newEdges);
		subvertices.set(newSubvertices);

		const initialKernel = initializeKernelState($vertices, $edges, alpha, beta);
		kernelData.set(initialKernel);
		previousEnergy.set(initialKernel.discreteEnergy || 0);

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
		cleanupInteractions = setupInteractions(
			$config.dim === 2 ? graphCanvas2D : graphCanvas3D,
			$vertices,
			updateVisualization
		);

		canvasTransform.set({ offsetX: 0, offsetY: 0, zoom: 1.0 });
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
		kernelData.set({ ...updatedKernel, disjointPairs: $kernelData.disjointPairs });
		updateVisualization();
	}

	function getEnergyChangeColor() {
		return $energyChange < 0 ? 'green' : 'red';
	}

	function switchDimension(event) {
		config.set({ ...$config, dim: parseInt(event.target.value) });
		cleanupInteractions();
		cleanupInteractions = setupInteractions(
			$config.dim === 2 ? graphCanvas2D : graphCanvas3D,
			$vertices,
			updateVisualization
		);
		regenerateGraph();
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
		<div>
			<label>
				Dimension:
				<select on:change={switchDimension}>
					<option value="2" selected>2D</option>
					<option value="3">3D</option>
				</select>
			</label>
		</div>
		<div class="energy-value">
			<p>Discrete Energy: {($discreteEnergy ?? 0).toFixed(4)}</p>
			<p style="color: {getEnergyChangeColor()}">
				Energy Change: {($energyChange ?? 0).toFixed(4)}
			</p>
		</div>
	</div>
	<div class="graph-section">
		<div class="graph-container" style="position: relative; width: {width}px; height: {height}px;">
			<canvas
				bind:this={graphCanvas2D}
				{width}
				{height}
				style="position: absolute; top: 0; left: 0; display: {$config.dim === 2
					? 'block'
					: 'none'};"
			></canvas>
			<canvas
				bind:this={graphCanvas3D}
				{width}
				{height}
				style="position: absolute; top: 0; left: 0; display: {$config.dim === 3
					? 'block'
					: 'none'};"
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
	.graph-container {
		position: relative;
	}
</style>
