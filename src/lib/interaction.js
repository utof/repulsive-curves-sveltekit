// src/lib/interaction.js
import { get, writable } from 'svelte/store';
import { canvasTransform } from '$lib/stores';

export function setupInteractions(canvas, vertices, updateFn) {
	let draggingVertex = null;
	let dragOffsetX = 0;
	let dragOffsetY = 0;
	let isDraggingCanvas = false; // For panning the canvas
	let lastMouseX = 0;
	let lastMouseY = 0;
	let isSpacePressed = false; // Track Space key state
	let isCtrlPressed = false; // Track Ctrl key state

	const MIN_ZOOM = 0.1; // Minimum zoom level
	const MAX_ZOOM = 5.0; // Maximum zoom level

	function getWorldCoords(screenX, screenY) {
		const { offsetX, offsetY, zoom } = get(canvasTransform);
		return {
			worldX: (screenX - offsetX) / zoom,
			worldY: (screenY - offsetY) / zoom
		};
	}

	function getScreenCoords(worldX, worldY) {
		const { offsetX, offsetY, zoom } = get(canvasTransform);
		return {
			screenX: worldX * zoom + offsetX,
			screenY: worldY * zoom + offsetY
		};
	}

	function handleMouseDown(event) {
		const rect = canvas.getBoundingClientRect();
		const mouseX = event.clientX - rect.left;
		const mouseY = event.clientY - rect.top;

		// Check if Space is pressed for canvas dragging (panning)
		if (isSpacePressed) {
			isDraggingCanvas = true;
			lastMouseX = mouseX;
			lastMouseY = mouseY;
			return; // Skip vertex dragging if panning
		}

		// Check for vertex dragging (existing behavior)
		const { zoom } = get(canvasTransform);
		for (let i = 0; i < vertices.length; i++) {
			const [vx, vy] = vertices[i];
			const { screenX, screenY } = getScreenCoords(vx, vy);
			const distance = Math.sqrt((mouseX - screenX) ** 2 + (mouseY - screenY) ** 2);
			if (distance <= 10 / zoom) {
				// Increased click radius for better hit detection
				draggingVertex = i;
				const worldCoords = getWorldCoords(mouseX, mouseY);
				dragOffsetX = vx - worldCoords.worldX;
				dragOffsetY = vy - worldCoords.worldY;
				break;
			}
		}
	}

	function handleMouseMove(event) {
		const rect = canvas.getBoundingClientRect();
		const mouseX = event.clientX - rect.left;
		const mouseY = event.clientY - rect.top;

		if (isDraggingCanvas && isSpacePressed) {
			// Panning: Update canvas offset with throttling for performance
			const dx = mouseX - lastMouseX;
			const dy = mouseY - lastMouseY;
			canvasTransform.update((transform) => ({
				...transform,
				offsetX: transform.offsetX + dx,
				offsetY: transform.offsetY + dy
			}));
			lastMouseX = mouseX;
			lastMouseY = mouseY;
			requestAnimationFrame(updateFn); // Use requestAnimationFrame for smoother updates
			return; // Skip vertex dragging if panning
		}

		if (draggingVertex !== null) {
			// Vertex dragging (optimized for performance)
			const worldCoords = getWorldCoords(mouseX, mouseY);
			let newX = worldCoords.worldX + dragOffsetX;
			let newY = worldCoords.worldY + dragOffsetY;

			// Keep vertices within canvas bounds (optional, adjust as needed)
			const { zoom } = get(canvasTransform);
			newX = Math.max(0, Math.min(canvas.width / zoom, newX));
			newY = Math.max(0, Math.min(canvas.height / zoom, newY));

			vertices[draggingVertex] = [newX, newY];
			requestAnimationFrame(updateFn); // Use requestAnimationFrame for smoother updates
		}
	}

	function handleMouseUp() {
		draggingVertex = null;
		isDraggingCanvas = false;
	}

	function handleMouseWheel(event) {
		if (isCtrlPressed) {
			event.preventDefault();
			const rect = canvas.getBoundingClientRect();
			const mouseX = event.clientX - rect.left;
			const mouseY = event.clientY - rect.top;

			// Calculate zoom factor
			const zoomFactor = event.deltaY < 0 ? 1.1 : 0.9; // Zoom in/out by 10%
			canvasTransform.update((transform) => {
				const newZoom = Math.min(MAX_ZOOM, Math.max(MIN_ZOOM, transform.zoom * zoomFactor));
				const worldCoords = getWorldCoords(mouseX, mouseY);
				const newOffsetX = mouseX - worldCoords.worldX * newZoom;
				const newOffsetY = mouseY - worldCoords.worldY * newZoom;
				return {
					offsetX: newOffsetX,
					offsetY: newOffsetY,
					zoom: newZoom
				};
			});

			requestAnimationFrame(updateFn); // Use requestAnimationFrame for smoother updates
		}
	}

	function handleKeyDown(event) {
		if (event.code === 'Space') {
			isSpacePressed = true;
		} else if (event.code === 'ControlLeft' || event.code === 'ControlRight') {
			isCtrlPressed = true;
		}
	}

	function handleKeyUp(event) {
		if (event.code === 'Space') {
			isSpacePressed = false;
			isDraggingCanvas = false; // Stop panning when Space is released
		} else if (event.code === 'ControlLeft' || event.code === 'ControlRight') {
			isCtrlPressed = false;
		}
	}

	// Add event listeners
	canvas.addEventListener('mousedown', handleMouseDown);
	canvas.addEventListener('mousemove', handleMouseMove);
	canvas.addEventListener('mouseup', handleMouseUp);
	canvas.addEventListener('mouseleave', handleMouseUp);
	canvas.addEventListener('wheel', handleMouseWheel, { passive: false });
	document.addEventListener('keydown', handleKeyDown);
	document.addEventListener('keyup', handleKeyUp);

	// Return a cleanup function
	return () => {
		canvas.removeEventListener('mousedown', handleMouseDown);
		canvas.removeEventListener('mousemove', handleMouseMove);
		canvas.removeEventListener('mouseup', handleMouseUp);
		canvas.removeEventListener('mouseleave', handleMouseUp);
		canvas.removeEventListener('wheel', handleMouseWheel);
		document.removeEventListener('keydown', handleKeyDown);
		document.removeEventListener('keyup', handleKeyUp);
	};
}
