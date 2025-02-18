<script>
	import { onMount } from 'svelte';
	import * as math from 'mathjs';

	let canvas;
	let ctx;
	let circles = [];
	let edges = [];
	let numberOfCircles = 0;
	const alpha = 2;
	const beta = 3;
	const learningRate = 0.9;
	let animationFrame;

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

	function energyGradient(circleIndex) {
		let grad = [0, 0];

		for (let [i, j] of edges) {
			if (i !== circleIndex && j !== circleIndex) continue;

			let y_i = [circles[i].x, circles[i].y];
			let y_j = [circles[j].x, circles[j].y];
			let diff = math.subtract(y_j, y_i);
			let distance = math.norm(diff);

			if (distance === 0) continue;

			let crossProduct = Math.abs(diff[0] * diff[1]);
			let k_value = Math.pow(crossProduct, alpha) / Math.pow(distance, beta);

			let gradFactor =
				(alpha * Math.pow(crossProduct, alpha - 1) * Math.sign(diff[0] * diff[1])) /
					Math.pow(distance, beta) -
				(beta * Math.pow(crossProduct, alpha)) / Math.pow(distance, beta + 1);

			let unitVector = math.divide(diff, distance);
			let gradContribution = math.multiply(unitVector, gradFactor);

			if (i === circleIndex) {
				grad = math.add(grad, gradContribution);
			} else {
				grad = math.subtract(grad, gradContribution);
			}
		}

		return grad;
	}

	function updateCirclePositions() {
		circles.forEach((circle, index) => {
			let grad = energyGradient(index);
			circle.x -= learningRate * grad[0];
			circle.y -= learningRate * grad[1];
		});
	}

	function computeDiscreteEnergy() {
		let energy = 0;
		for (let [i, j] of edges) {
			let y_i = [circles[i].x, circles[i].y];
			let y_j = [circles[j].x, circles[j].y];
			let segmentLength = math.distance(y_i, y_j);
			let k_value = k_alpha_beta(y_i, y_j, alpha, beta);
			energy += k_value * segmentLength * segmentLength;
		}
		return energy;
	}

	function draw() {
		if (!ctx) return;
		ctx.clearRect(0, 0, canvas.width, canvas.height);

		// Draw edges
		ctx.strokeStyle = '#ccc';
		edges.forEach(([i, j]) => {
			ctx.beginPath();
			ctx.moveTo(circles[i].x, circles[i].y);
			ctx.lineTo(circles[j].x, circles[j].y);
			ctx.stroke();
		});

		// Draw circles
		circles.forEach((circle) => {
			ctx.beginPath();
			ctx.arc(circle.x, circle.y, circle.radius, 0, 2 * Math.PI);
			ctx.fillStyle = circle.color;
			ctx.fill();
		});
	}

	function animate() {
		updateCirclePositions();
		draw();
		animationFrame = requestAnimationFrame(animate);
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
		draw();
		animate();

		window.addEventListener('resize', () => {
			resizeCanvas();
			generateCircles();
			draw();
		});
	});
</script>

<svelte:head>
	<title>Gradient-Based Energy Minimization</title>
</svelte:head>

<canvas id="myCanvas"></canvas>

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
