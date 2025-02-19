<script>
	import { onMount, onDestroy } from 'svelte';
	import * as math from 'mathjs';
	import { getColor, getRadius } from '$lib/graphUtils';
	import { discreteTangentPointEnergy, calculateEdgeProperties } from '$lib/energyCalculations';

	let canvas;
	let ctx;
	let vertices = [];
	let edges = [];
	let energy = 0;
	let vertexData = [];
	let edgeProps = { edgeLengths: [], edgeTangents: [], edgeMidpoints: [] };
	let animationFrameId;
	let isAnimating = false;

	const width = 1920;
	const height = 1080;
	const alpha = 3;
	const beta = 6;
	const learningRate = 0.01; // Adjust for animation speed
	const maxIterations = 100; // Maximum iterations per step
	const TANGENT_ARROW_LENGTH = 20; // Length of tangent arrow visualization
	let currentIteration = 0;

	onMount(() => {
		canvas = document.getElementById('graphCanvas');
		ctx = canvas.getContext('2d');
		generateRandomGraph();
		calculateEnergy();
		drawGraph();
	});

	onDestroy(() => {
		if (animationFrameId) {
			// Check if animationFrameId is defined
			cancelAnimationFrame(animationFrameId); // Clean up animation frame
		}
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
				if (Math.random() < 0.2) {
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
		currentIteration = 0; // Reset iteration count on regeneration
	}

	function calculateEnergy() {
		const result = discreteTangentPointEnergy(vertices, edges, alpha, beta);
		energy = result.totalEnergy;
		vertexData = result.vertexData;
		edgeProps = calculateEdgeProperties(vertices, edges);
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

		// Draw edges
		ctx.strokeStyle = 'gray';
		ctx.lineWidth = 1;
		for (const edge of edges) {
			ctx.beginPath();
			ctx.moveTo(vertices[edge[0]][0], vertices[edge[0]][1]);
			ctx.lineTo(vertices[edge[1]][0], vertices[edge[1]][1]);
			ctx.stroke();
		}

		// Draw vertices with numbers
		for (let i = 0; i < vertices.length; i++) {
			const vertexEnergy = vertexData[i]?.energy || 0;
			const radius = getRadius(vertexEnergy, vertexData);

			// Draw vertex circle
			ctx.beginPath();
			ctx.arc(vertices[i][0], vertices[i][1], radius, 0, 2 * Math.PI);
			ctx.fillStyle = getColor(vertexEnergy, vertexData);
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

	function animate() {
		if (!isAnimating) return;

		if (currentIteration >= maxIterations) {
			isAnimating = false; // Stop after maxIterations
			return;
		}
		currentIteration++;

		calculateEnergy(); // Recalculate energy and gradients

		// Update vertex positions based on gradients
		for (let i = 0; i < vertices.length; i++) {
			vertices[i] = math.subtract(vertices[i], math.multiply(learningRate, vertexData[i].gradient));

			// Keep vertices within bounds
			vertices[i][0] = Math.max(0, Math.min(width, vertices[i][0]));
			vertices[i][1] = Math.max(0, Math.min(height, vertices[i][1]));
		}

		drawGraph();
		animationFrameId = requestAnimationFrame(animate);
	}

	function startAnimation() {
		if (!isAnimating) {
			isAnimating = true;
			currentIteration = 0; // Reset iteration count
			animate();
		}
	}

	function stopAnimation() {
		isAnimating = false;
	}

	function regenerateAndAnimate() {
		generateRandomGraph();
		startAnimation();
	}
</script>

<canvas id="graphCanvas" {width} {height}></canvas>

<p>Total Energy: {energy.toFixed(2)}</p>
<div class="vertex-data">
	<h3>Vertex Data:</h3>
	{#each vertexData as data, i}
		<div class="vertex-info">
			<p>Vertex {i}:</p>
			<ul>
				<li>Energy: {data.energy.toFixed(4)}</li>
				<li>Gradient: [{data.gradient.map((g) => g.toFixed(4)).join(', ')}]</li>
				<li>Count: {data.count}</li>
			</ul>
		</div>
	{/each}
</div>

<button on:click={regenerateAndAnimate}>Regenerate and Animate</button>
<button on:click={startAnimation}>Start Animation</button>
<button on:click={stopAnimation}>Stop Animation</button>

<style>
	.vertex-data {
		margin-top: 20px;
		padding: 10px;
		border: 1px solid #ccc;
		border-radius: 4px;
		max-height: 300px;
		overflow-y: auto;
	}

	.vertex-info {
		margin: 10px 0;
		padding: 5px;
		border-bottom: 1px solid #eee;
	}

	.vertex-info p {
		margin: 0;
		font-weight: bold;
	}

	.vertex-info ul {
		margin: 5px 0;
		padding-left: 20px;
	}

	.vertex-info li {
		margin: 2px 0;
	}
</style>
