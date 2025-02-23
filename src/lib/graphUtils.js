/**
 * Calculate color based on vertex energy
 * @param {number} vertexEnergy - Energy value for the vertex
 * @param {Array} vertexData - Array of vertex data containing energy values
 * @returns {string} RGB color string
 */
export function getColor(vertexEnergy, vertexData) {
	const maxEnergy = Math.max(...vertexData.map((v) => v.energy), 1e-9); // Avoid division by zero
	const normalizedEnergy = vertexEnergy / maxEnergy;
	const r = Math.floor(255 * normalizedEnergy);
	const g = 0;
	const b = Math.floor(255 * (1 - normalizedEnergy));
	return `rgb(${r}, ${g}, ${b})`;
}

/**
 * Calculate radius based on vertex energy
 * @param {number} vertexEnergy - Energy value for the vertex
 * @param {Array} vertexData - Array of vertex data containing energy values
 * @returns {number} Radius value
 */
export function getRadius(vertexEnergy, vertexData) {
	const maxEnergy = Math.max(...vertexData.map((v) => v.energy), 1e-9); // Avoid division by zero
	const normalizedEnergy = vertexEnergy / maxEnergy;
	return 2 + 8 * normalizedEnergy;
}

export function drawArrow(ctx, fromX, fromY, dirX, dirY, length) {
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

export function generateRandomGraph(width, height) {
	const numVertices = Math.floor(Math.random() * 10) + 5;
	const vertices = [];
	const edges = [];

	for (let i = 0; i < numVertices; i++) {
		vertices.push([
			Math.random() * width, // x within canvas width
			Math.random() * height // y within canvas height
		]);
	}

	const edgeSet = new Set();
	for (let i = 0; i < numVertices; i++) {
		for (let j = i + 1; j < numVertices; j++) {
			if (Math.random() < 0.1) {
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
	return { vertices, edges };
}
