import * as math from 'mathjs';

// Function to compute k_{β}^{α}(yi, yj, TI)
function k_alpha_beta(y_i, y_j, T_I, alpha, beta) {
	const diff = math.subtract(y_j, y_i);
	const crossProduct = math.norm(math.cross(T_I, diff));
	const distance = math.norm(diff);
	return Math.pow(crossProduct, alpha) / Math.pow(distance, beta);
}

// Compute the discrete tangent point energy
function computeDiscreteEnergy(edges, points, tangents, segmentLengths, alpha, beta) {
	let energy = 0;

	for (let I of edges) {
		for (let J of edges) {
			if (I !== J && I.filter((i) => J.includes(i)).length === 0) {
				// Ensure disjoint segments
				let sum_k = 0;

				for (let i of I) {
					for (let j of J) {
						sum_k += k_alpha_beta(points[i], points[j], tangents[I], alpha, beta);
					}
				}

				energy += (1 / 4) * sum_k * segmentLengths[I] * segmentLengths[J];
			}
		}
	}

	return energy;
}

// Example usage
const edges = [
	[0, 1],
	[2, 3]
]; // Example edge indices
const points = [
	[0, 0, 0],
	[1, 0, 0],
	[2, 1, 0],
	[3, 1, 0]
]; // Example 3D points
tangents = [
	[0, 1, 0],
	[0, 1, 0]
]; // Example tangents
const segmentLengths = { 0: 1, 1: 1 }; // Example segment lengths
const alpha = 2;
const beta = 3;

const energy = computeDiscreteEnergy(edges, points, tangents, segmentLengths, alpha, beta);
console.log('Discrete Tangent Point Energy:', energy);
