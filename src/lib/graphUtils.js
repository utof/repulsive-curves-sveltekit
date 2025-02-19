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
