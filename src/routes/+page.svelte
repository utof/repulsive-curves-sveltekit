<!-- +page.svelte -->
<script>
	import { onMount } from 'svelte';
	import * as math from 'mathjs';

	let vertices = [];
	let edges = [];
	let energy = 0;
	let vertexData = []; // Store energy-related data per vertex

	onMount(() => {
		generateRandomGraph();
		calculateEnergy();
	});

	function generateRandomGraph() {
		const numVertices = Math.floor(Math.random() * (50 - 5 + 1)) + 5; // 5 to 50 vertices
		vertices = [];
		edges = [];

		// Generate vertices with random 3D coordinates
		for (let i = 0; i < numVertices; i++) {
			vertices.push([
				Math.random() * 2 - 1, // x: -1 to 1
				Math.random() * 2 - 1, // y: -1 to 1
				Math.random() * 2 - 1 // z: -1 to 1
			]);
		}

		// Generate edges (at most one edge per vertex pair)
		const edgeSet = new Set(); // Use a Set to prevent duplicate edges
		for (let i = 0; i < numVertices; i++) {
			for (let j = i + 1; j < numVertices; j++) {
				if (Math.random() < 0.05) {
					// 30% chance of creating an edge
					const edge = [i, j];
					const edgeString = `${i}-${j}`; // Create a unique string for the edge
					if (!edgeSet.has(edgeString)) {
						edges.push(edge);
						edgeSet.add(edgeString);
					}
				}
			}
		}
	}

	function discreteTangentPointEnergy(vertices, edges, alpha, beta) {
		let totalEnergy = 0;
		vertexData = vertices.map(() => ({ energy: 0, count: 0 })); // Initialize per-vertex data

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
				const tx_cross_p_minus_q = math.cross(ti, p_minus_q);

				const numerator = Math.pow(math.norm(tx_cross_p_minus_q), beta);
				const denominator = Math.pow(math.norm(p_minus_q), beta);

				const k_beta = numerator / denominator;
				const edgeEnergy = k_beta * edgeLengths[i] * edgeLengths[j];
				totalEnergy += edgeEnergy;

				// Accumulate energy contributions per vertex
				vertexData[edges[i][0]].energy += edgeEnergy;
				vertexData[edges[i][0]].count += 1;
				vertexData[edges[i][1]].energy += edgeEnergy;
				vertexData[edges[i][1]].count += 1;
			}
		}
		// Average the energy contributions for each vertex
		for (let i = 0; i < vertexData.length; i++) {
			if (vertexData[i].count > 0) {
				vertexData[i].energy /= vertexData[i].count;
			}
		}

		return totalEnergy;
	}

	function calculateEnergy() {
		const alpha = 3;
		const beta = 6;
		energy = discreteTangentPointEnergy(vertices, edges, alpha, beta);
	}

	function regenerateGraph() {
		generateRandomGraph();
		calculateEnergy();
	}

	function getColor(vertexEnergy) {
		// Normalize energy to [0, 1] range (adjust maxEnergy as needed)
		const maxEnergy = Math.max(...vertexData.map((v) => v.energy)); // Find maximum vertex energy
		const normalizedEnergy = maxEnergy > 0 ? vertexEnergy / maxEnergy : 0;

		// Simple color mapping: blue (low energy) to red (high energy)
		const r = Math.floor(255 * normalizedEnergy);
		const g = 0;
		const b = Math.floor(255 * (1 - normalizedEnergy));
		return `rgb(${r}, ${g}, ${b})`;
	}

	function getRadius(vertexEnergy) {
		// Normalize energy to [0, 1] range (adjust maxEnergy as needed)
		const maxEnergy = Math.max(...vertexData.map((v) => v.energy)); // Find maximum vertex energy
		const normalizedEnergy = maxEnergy > 0 ? vertexEnergy / maxEnergy : 0;
		// Map normalized energy to a radius between 2 and 10 (adjust as needed)
		return 2 + 8 * normalizedEnergy;
	}
</script>

<div class="container">
	{#each vertices as vertex, i}
		<div
			class="vertex"
			style="
          left: {vertex[0] * 200 + 250}px; /* Scale and center */
          top: {vertex[1] * 200 + 250}px;
          background-color: {getColor(vertexData[i]?.energy || 0)};
          width: {getRadius(vertexData[i]?.energy || 0)}px;
          height: {getRadius(vertexData[i]?.energy || 0)}px;
        "
		>
			<!-- Optional: Display vertex index or other info -->
			<!-- {i} -->
		</div>
	{/each}

	{#each edges as edge}
		{@const v1 = vertices[edge[0]]}
		{@const v2 = vertices[edge[1]]}
		{@const dx = v2[0] - v1[0]}
		{@const dy = v2[1] - v1[1]}
		{@const length = Math.sqrt(dx * dx + dy * dy) * 200}
		{@const angle = (Math.atan2(dy, dx) * 180) / Math.PI}
		<div
			class="edge"
			style="
          left: {v1[0] * 200 + 250}px;
          top: {v1[1] * 200 + 250}px;
          width: {length}px;
          transform: rotate({angle}deg);
        "
		></div>
	{/each}
</div>

<p>Total Energy: {energy.toFixed(2)}</p>
<button on:click={regenerateGraph}>Regenerate Graph</button>

<style>
	.container {
		width: 100%;
		height: 500px; /* Adjust as needed */
		position: relative;
		border: 1px solid black;
	}
	.vertex {
		position: absolute;
		border-radius: 50%; /* Make it a circle */
		transform: translate(-50%, -50%); /* Center the circle */
	}
	.edge {
		position: absolute;
		background-color: gray;
		height: 1px; /* Thin line */
		transform-origin: 0 50%; /* Rotate around the starting point */
	}
</style>
