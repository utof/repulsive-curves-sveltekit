<script>
	import { onMount } from 'svelte';

	let canvas;
	let ctx;
	let circles = [];
	let numberOfCircles = 0;

	function resizeCanvas() {
		canvas.width = window.innerWidth;
		canvas.height = window.innerHeight;
	}

	function generateCircles() {
		numberOfCircles = Math.floor(Math.random() * 4) + 2; // Random number of circles between 2 and 5
		circles = [];
		for (let i = 0; i < numberOfCircles; i++) {
			circles.push({
				x: Math.random() * canvas.width,
				y: Math.random() * canvas.height,
				radius: Math.random() * 20 + 10, // Random radius between 10 and 30
				color: `hsl(${Math.random() * 360}, 70%, 50%)` // Random color
			});
		}
	}

	function draw() {
		if (!ctx) return;
		ctx.clearRect(0, 0, canvas.width, canvas.height);

		// Draw lines
		for (let i = 0; i < circles.length; i++) {
			for (let j = i + 1; j < circles.length; j++) {
				if (Math.random() < 0.7) {
					// Randomly draw lines (70% probability)
					ctx.beginPath();
					ctx.moveTo(circles[i].x, circles[i].y);
					ctx.lineTo(circles[j].x, circles[j].y);
					ctx.strokeStyle = '#ccc';
					ctx.stroke();
				}
			}
		}

		// Draw circles
		circles.forEach((circle) => {
			ctx.beginPath();
			ctx.arc(circle.x, circle.y, circle.radius, 0, 2 * Math.PI);
			ctx.fillStyle = circle.color;
			ctx.fill();
		});
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

		// Redraw on window resize
		window.addEventListener('resize', () => {
			resizeCanvas();
			generateCircles();
			draw();
		});
	});
</script>

<svelte:head>
	<title>Random Circles</title>
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
