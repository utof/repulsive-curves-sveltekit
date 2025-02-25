// src/lib/graphUtils.js
import { get } from 'svelte/store';
import { config } from './stores';

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
    ctx.lineTo(toX - headLength * Math.cos(angle - headAngle), toY - headLength * Math.sin(angle - headAngle));
    ctx.moveTo(toX, toY);
    ctx.lineTo(toX - headLength * Math.cos(angle + headAngle), toY - headLength * Math.sin(angle + headAngle));
    ctx.stroke();
}

export function generateSubvertices(vertices, edges) {
    const dim = get(config).dim;
    const subvertices = [];
    const { subvertexGap } = get(config);

    for (const edge of edges) {
        const [v1, v2] = edge;
        const p1 = vertices[v1];
        const p2 = vertices[v2];
        const diff = p2.map((coord, i) => coord - p1[i]);
        const edgeLength = Math.sqrt(diff.reduce((sum, d) => sum + d * d, 0));
        const subvertexCount = Math.max(1, Math.floor(edgeLength / subvertexGap));

        for (let i = 1; i <= subvertexCount; i++) {
            const t = i / (subvertexCount + 1);
            const pos = p1.map((coord, idx) => coord + t * diff[idx]);
            subvertices.push({ position: pos, edge });
        }
    }
    // return subvertices;
    return [];
}

export function generateRandomGraph(width, height) {
    const dim = get(config).dim;
    const numVertices = Math.floor(Math.random() * 10) + 5;
    const vertices = [];

    for (let i = 0; i < numVertices; i++) {
        const vertex = [];
        for (let d = 0; d < dim; d++) {
            vertex.push(Math.random() * (d < 2 ? width : height));
        }
        vertices.push(vertex);
    }

    const edges = [];
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
    const dim = get(config).dim;
    const vertices = [];
    const edges = [];

    for (let i = 0; i < 3; i++) {
        const vertex = new Array(dim).fill(0);
        vertex[0] = width * 0.3;
        vertex[1] = height * (0.2 + i * 0.3);
        if (dim === 3) vertex[2] = Math.random() * height;
        vertices.push(vertex);
    }

    for (let i = 0; i < 3; i++) {
        const vertex = new Array(dim).fill(0);
        vertex[0] = width * 0.7;
        vertex[1] = height * (0.2 + i * 0.3);
        if (dim === 3) vertex[2] = Math.random() * height;
        vertices.push(vertex);
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
    const dim = get(config).dim;
    const vertices = [];
    const edges = [];

    for (let i = 0; i < 2; i++) {
        const vertex = new Array(dim).fill(0);
        vertex[0] = width * 0.3;
        vertex[1] = height * (0.3 + i * 0.4);
        if (dim === 3) vertex[2] = Math.random() * height;
        vertices.push(vertex);
    }

    for (let i = 0; i < 3; i++) {
        const vertex = new Array(dim).fill(0);
        vertex[0] = width * 0.7;
        vertex[1] = height * (0.2 + i * 0.3);
        if (dim === 3) vertex[2] = Math.random() * height;
        vertices.push(vertex);
    }

    for (let i = 0; i < 2; i++) {
        for (let j = 0; j < 3; j++) {
            edges.push([i, 2 + j]);
        }
    }

    const subvertices = generateSubvertices(vertices, edges);
    return { vertices, edges, subvertices };
}