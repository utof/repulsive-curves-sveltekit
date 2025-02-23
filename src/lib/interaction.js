// src/lib/interaction.js

export function setupInteractions(canvas, vertices, updateFn) {
	let draggingVertex = null;
	let dragOffsetX = 0;
	let dragOffsetY = 0;

	function handleMouseDown(event) {
		const rect = canvas.getBoundingClientRect();
		const mouseX = event.clientX - rect.left;
		const mouseY = event.clientY - rect.top;

		for (let i = 0; i < vertices.length; i++) {
			const [vx, vy] = vertices[i];
			const distance = Math.sqrt((mouseX - vx) ** 2 + (mouseY - vy) ** 2);
			if (distance <= 5) {
				draggingVertex = i;
				dragOffsetX = vx - mouseX;
				dragOffsetY = vy - mouseY;
				break;
			}
		}
	}

	function handleMouseMove(event) {
		if (draggingVertex !== null) {
			const rect = canvas.getBoundingClientRect();
			const mouseX = event.clientX - rect.left;
			const mouseY = event.clientY - rect.top;

			let newX = mouseX + dragOffsetX;
			let newY = mouseY + dragOffsetY;

			vertices[draggingVertex] = [newX, newY];
			updateFn(); // Call the update function
		}
	}

	function handleMouseUp() {
		draggingVertex = null;
	}

	canvas.addEventListener('mousedown', handleMouseDown);
	canvas.addEventListener('mousemove', handleMouseMove);
	canvas.addEventListener('mouseup', handleMouseUp);
	canvas.addEventListener('mouseleave', handleMouseUp);

	// Return a cleanup function
	return () => {
		canvas.removeEventListener('mousedown', handleMouseDown);
		canvas.removeEventListener('mousemove', handleMouseMove);
		canvas.removeEventListener('mouseup', handleMouseUp);
		canvas.removeEventListener('mouseleave', handleMouseUp);
	};
}
