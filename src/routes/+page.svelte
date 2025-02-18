<script>
	import { onMount } from 'svelte';
	import * as math from 'mathjs';

	let canvas;
	let ctx;
	let circles = [];
	let edges = [];
	let numberOfCircles = 0;
	let energy = 0;
	const alpha = 2;
	const beta = 3;

	function resizeCanvas() {
		canvas.width = window.innerWidth;
		canvas.height = window.innerHeight;
	}

	function generateCircles() {
		numberOfCircles = Math.floor(Math.random() * 4) + 2;
		circles = [];
		edges = [];
		for (let i = 0; i < numberOfCircles; i++) {
			circles.push({
				x: Math.random() * canvas.width,
				y: Math.random() * canvas.height,
				radius: Math.random() * 20 + 10,
				color: `hsl(${Math.random() * 360}, 70%, 50%)`
			});
		}

		for (let i = 0; i < circles.length; i++) {
			for (let j = i + 1; j < circles.length; j++) {
				if (Math.random() < 0.7) {
					edges.push([i, j]);
				}
			}
		}
	}

	function k_alpha_beta(y_i, y_j, alpha, beta) {
		const diff = math.subtract(y_j, y_i);
		const crossProduct = Math.abs(diff[0] * diff[1]);
		const distance = math.norm(diff);
		return Math.pow(crossProduct, alpha) / Math.pow(distance, beta);
	}

	function computeDiscreteEnergy() {
		energy = 0;
		console.log('Computing energy...');

		for (let [i, j] of edges) {
			let y_i = [circles[i].x, circles[i].y];
			let y_j = [circles[j].x, circles[j].y];
			let segmentLength = math.distance(y_i, y_j);
			let k_value = k_alpha_beta(y_i, y_j, alpha, beta);
			energy += k_value * segmentLength * segmentLength;
			console.log(
				`Edge (${i}, ${j}) -> Energy contribution:`,
				k_value * segmentLength * segmentLength
			);
		}

		console.log('Total Discrete Tangent Point Energy:', energy);
		visualizeEnergy();
	}

	function draw() {
		if (!ctx) return;
		ctx.clearRect(0, 0, canvas.width, canvas.height);

		ctx.strokeStyle = '#ccc';
		edges.forEach(([i, j]) => {
			ctx.beginPath();
			ctx.moveTo(circles[i].x, circles[i].y);
			ctx.lineTo(circles[j].x, circles[j].y);
			ctx.stroke();
		});

		circles.forEach((circle) => {
			ctx.beginPath();
			ctx.arc(circle.x, circle.y, circle.radius, 0, 2 * Math.PI);
			ctx.fillStyle = circle.color;
			ctx.fill();
		});
	}

	function visualizeEnergy() {
		const energyCanvas = document.getElementById('energyCanvas');
		const ctxEnergy = energyCanvas.getContext('2d');

		energyCanvas.width = window.innerWidth;
		energyCanvas.height = window.innerHeight;

		ctxEnergy.clearRect(0, 0, energyCanvas.width, energyCanvas.height);

		const maxRadius = 100;
		const colorIntensity = Math.min(255, Math.floor(energy * 10));

		ctxEnergy.beginPath();
		ctxEnergy.arc(
			energyCanvas.width / 2,
			energyCanvas.height / 2,
			Math.min(maxRadius, energy * 2),
			0,
			2 * Math.PI
		);
		ctxEnergy.fillStyle = `rgb(${colorIntensity}, 50, 150)`;
		ctxEnergy.fill();
	}

	onMount(() => {
		canvas = document.getElementById('myCanvas');
		ctx = canvas.getContext('2d');
		if (!ctx) {
			console.error('Canvas context not available');
			return;
		}
		resizeCanvas();
		generateCircles();
		computeDiscreteEnergy();
		draw();

		window.addEventListener('resize', () => {
			resizeCanvas();
			generateCircles();
			computeDiscreteEnergy();
			draw();
		});
	});
</script>

<svelte:head>
	<title>Random Circles Energy</title>
</svelte:head>

<canvas id="myCanvas"></canvas>
<canvas
	id="energyCanvas"
	style="position: fixed; top: 0; left: 0; width: 100%; height: 100%; pointer-events: none;"
></canvas>

<style>
	canvas {
		position: fixed;
		top: 0;
		left: 0;
		width: 100%;
		height: 100%;
		z-index: -1;
	}
</style>
