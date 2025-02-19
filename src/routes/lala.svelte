<!-- +page.svelte -->
<script>
	import { onMount, onDestroy } from 'svelte';
	import * as math from 'mathjs';

	let canvas;
	let ctx;
	let vertices = [];
	let edges = [];
	let energy = 0;
	let vertexData = [];
	let animationFrameId;
	let isAnimating = false;

	const width = 500;
	const height = 500;
	const alpha = 3;
	const beta = 6;
	const learningRate = 0.01; // Adjust for animation speed
	const maxIterations = 100; // Maximum iterations per step
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

	function discreteTangentPointEnergy(vertices, edges, alpha, beta) {
		let totalEnergy = 0;
		vertexData = vertices.map(() => ({ energy: 0, count: 0, gradient: [0, 0] }));

		const edgeLengths = edges.map((edge) => math.distance(vertices[edge[0]], vertices[edge[1]]));
		const edgeTangents = edges.map((edge, i) =>
			math.divide(math.subtract(vertices[edge[1]], vertices[edge[0]]), edgeLengths[i])
		);

		for (let i = 0; i < edges.length; i++) {
			for (let j = 0; j < edges.length; j++) {
				if (i === j || edges[i].some((v) => edges[j].includes(v))) continue;

				const ti = edgeTangents[i];
				const xi = vertices[edges[i][0]];
				const xj = vertices[edges[j][0]];

				const p_minus_q = math.subtract(xi, xj);
				const tx_cross_p_minus_q = math.cross([...ti, 0], [...p_minus_q, 0]); // 2D cross product (z-component)

				const numerator = Math.pow(Math.abs(tx_cross_p_minus_q[2]), beta); // Use abs for 2D
				const denominator = Math.pow(math.norm(p_minus_q), beta);

				const k_beta = numerator / denominator;
				const edgeEnergy = k_beta * edgeLengths[i] * edgeLengths[j];
				totalEnergy += edgeEnergy;

				// --- Gradient Calculation (for animation) ---
				const d_k_beta_d_xi = math.multiply(
					(beta * k_beta) / math.norm(p_minus_q),
					math.subtract(
						math.multiply(math.dot(ti, p_minus_q) / math.norm(p_minus_q), ti),
						math.multiply(
							(tx_cross_p_minus_q[2] * tx_cross_p_minus_q[2]) / (numerator + 1e-9),
							p_minus_q
						) //add small number to prevent NaN
					)
				);

				// Accumulate gradient and energy per vertex
				vertexData[edges[i][0]].gradient = math.add(
					vertexData[edges[i][0]].gradient,
					math.multiply(edgeLengths[i] * edgeLengths[j], d_k_beta_d_xi)
				);
				vertexData[edges[i][1]].gradient = math.subtract(
					vertexData[edges[i][1]].gradient,
					math.multiply(edgeLengths[i] * edgeLengths[j], d_k_beta_d_xi)
				);

				vertexData[edges[i][0]].energy += edgeEnergy;
				vertexData[edges[i][0]].count += 1;
				vertexData[edges[i][1]].energy += edgeEnergy;
				vertexData[edges[i][1]].count += 1;
			}
		}

		for (let i = 0; i < vertexData.length; i++) {
			if (vertexData[i].count > 0) {
				vertexData[i].energy /= vertexData[i].count;
			}
		}

		return totalEnergy;
	}

	function calculateEnergy() {
		energy = discreteTangentPointEnergy(vertices, edges, alpha, beta);
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

		// Draw vertices
		for (let i = 0; i < vertices.length; i++) {
			const vertexEnergy = vertexData[i]?.energy || 0;
			ctx.beginPath();
			ctx.arc(vertices[i][0], vertices[i][1], getRadius(vertexEnergy), 0, 2 * Math.PI);
			ctx.fillStyle = getColor(vertexEnergy);
			ctx.fill();
		}
	}

	function getColor(vertexEnergy) {
		const maxEnergy = Math.max(...vertexData.map((v) => v.energy), 1e-9); // Avoid division by zero
		const normalizedEnergy = vertexEnergy / maxEnergy;
		const r = Math.floor(255 * normalizedEnergy);
		const g = 0;
		const b = Math.floor(255 * (1 - normalizedEnergy));
		return `rgb(${r}, ${g}, ${b})`;
	}

	function getRadius(vertexEnergy) {
		const maxEnergy = Math.max(...vertexData.map((v) => v.energy), 1e-9); // Avoid division by zero
		const normalizedEnergy = vertexEnergy / maxEnergy;
		return 2 + 8 * normalizedEnergy;
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
		regenerateGraph();
		startAnimation();
	}
</script>

<canvas id="graphCanvas" {width} {height}></canvas>

<p>Total Energy: {energy.toFixed(2)}</p>
<button on:click={regenerateAndAnimate}>Regenerate and Animate</button>
<button on:click={startAnimation}>Start Animation</button>
<button on:click={stopAnimation}>Stop Animation</button>
