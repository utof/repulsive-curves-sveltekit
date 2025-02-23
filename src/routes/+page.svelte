<script>
	import { onMount, onDestroy } from 'svelte';
	import {
		calculateEdgeProperties,
		calculateDisjointEdgePairs,
		calculateDiscreteKernel,
		calculateDiscreteEnergy
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
	let kernelMatrix = null;
	let disjointPairs = [];

	const width = 700;
	const height = 700;
	const alpha = 3;
	const beta = 6;
	const stepSize = 0.01;
	const maxIterations = 100;
	let optimizer = null;

	onMount(() => {
		canvas = document.getElementById('graphCanvas');
		ctx = canvas.getContext('2d');
		kernelCanvas = document.getElementById('kernelCanvas');
		kernelCtx = kernelCanvas.getContext('2d');

		regenerateGraph(); // Initial graph

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
	});

	function calculateAndDrawKernel() {
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
		drawKernelMatrix(kernelCtx, kernelMatrix);
	}

	function regenerateGraph() {
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
