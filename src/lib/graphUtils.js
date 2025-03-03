// src/lib/graphUtils.js
import { get } from 'svelte/store';
import { config } from './stores';

export function getColor(vertexEnergy, vertexData) {
	const maxEnergy = Math.max(...vertexData.map((v) => v.energy), 1e-9);
	const normalizedEnergy = vertexEnergy / maxEnergy;
	const r = Math.floor(255 * normalizedEnergy);
	const g = 0;
	const b = Math.floor(255 * (1 - normalizedEnergy));
	return `rgb(${r}, ${g}, ${b})`;
}

export function getRadius(vertexEnergy, vertexData) {
	const maxEnergy = Math.max(...vertexData.map((v) => v.energy), 1e-9);
	const normalizedEnergy = vertexEnergy / maxEnergy;
	return 2 + 8 * normalizedEnergy;
}

export function drawArrow(ctx, fromX, fromY, dirX, dirY, length) {
	const headLength = 7;
	const headAngle = Math.PI / 6;

	const toX = fromX + dirX * length;
	const toY = fromY + dirY * length;

	ctx.beginPath();
	ctx.moveTo(fromX, fromY);
	ctx.lineTo(toX, toY);
	ctx.stroke();

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

export function generateSubvertices(vertices, edges) {
	const subvertices = [];
	const { subvertexGap } = get(config);

	for (const edge of edges) {
		const [v1, v2] = edge;
		const p1 = vertices[v1];
		const p2 = vertices[v2];
		const dx = p2[0] - p1[0];
		const dy = p2[1] - p1[1];
		const edgeLength = Math.sqrt(dx * dx + dy * dy);
		const subvertexCount = Math.max(1, Math.floor(edgeLength / subvertexGap));

		for (let i = 1; i <= subvertexCount; i++) {
			const t = i / (subvertexCount + 1);
			const x = p1[0] + t * dx;
			const y = p1[1] + t * dy;
			subvertices.push({ position: [x, y], edge });
		}
	}
	return subvertices;
}

export function generateRandomGraph(width, height) {
	const numVertices = Math.floor(Math.random() * 10) + 5;
	const vertices = [];
	const edges = [];

	for (let i = 0; i < numVertices; i++) {
		vertices.push([Math.random() * width, Math.random() * height]);
	}

	const edgeSet = new Set();
	for (let i = 0; i < numVertices; i++) {
		for (let j = i + 1; j < numVertices; j++) {
			if (Math.random() < 0.1) {
				const edge = [i, j];
				const edgeString = `${i}-${j}`;
				if (!edgeSet.has(edgeString)) {
					edges.push(edge);
					edgeSet.add(edgeString);
				}
			}
		}
	}

	const subvertices = generateSubvertices(vertices, edges);
	return { vertices, edges, subvertices };
}

export function generateBipartiteGraph(width, height) {
	const vertices = [];
	const edges = [];

	for (let i = 0; i < 3; i++) {
		vertices.push([width * 0.3, height * (0.2 + i * 0.3)]);
	}

	for (let i = 0; i < 3; i++) {
		vertices.push([width * 0.7, height * (0.2 + i * 0.3)]);
	}

	for (let i = 0; i < 3; i++) {
		for (let j = 0; j < 3; j++) {
			edges.push([i, 3 + j]);
		}
	}

	const subvertices = generateSubvertices(vertices, edges);
	return { vertices, edges, subvertices };
}

export function generate2x3BipartiteGraph(width, height) {
	const vertices = [];
	const edges = [];

	for (let i = 0; i < 2; i++) {
		vertices.push([width * 0.3, height * (0.3 + i * 0.4)]);
	}

	for (let i = 0; i < 3; i++) {
		vertices.push([width * 0.7, height * (0.2 + i * 0.3)]);
	}

	for (let i = 0; i < 2; i++) {
		for (let j = 0; j < 3; j++) {
			edges.push([i, 2 + j]);
		}
	}

	const subvertices = generateSubvertices(vertices, edges);
	return { vertices, edges, subvertices };
}

export function generateSimpleSquareGraph(width, height) {
	const vertices = [
		[width * 0.3, height * 0.3],
		[width * 0.7, height * 0.3],
		[width * 0.7, height * 0.7],
		[width * 0.3, height * 0.7]
	];

	const edges = [
		[0, 1],
		[1, 2],
		[2, 3],
		[3, 0]
	];

	const subvertices = generateSubvertices(vertices, edges);
	return { vertices, edges, subvertices };
}

export function generateCompleteSquareGraph(width, height) {
	const vertices = [
		[width * 0.3, height * 0.3],
		[width * 0.7, height * 0.3],
		[width * 0.7, height * 0.7],
		[width * 0.3, height * 0.7]
	];

	const edges = [
		[0, 1],
		[0, 2],
		[0, 3],
		[1, 2],
		[1, 3],
		[2, 3]
	];

	const subvertices = generateSubvertices(vertices, edges);
	return { vertices, edges, subvertices };
}