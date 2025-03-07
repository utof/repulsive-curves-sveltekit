// src/lib/graphUtils.js
import { get } from 'svelte/store';
import { config } from './stores';

/**
 * Projects a 3D point onto a 2D canvas
 * Simple orthographic projection
 * @param {Array} point - 3D point [x,y,z]
 * @returns {Array} - 2D point [x,y]
 */
export function project3Dto2D(point) {
    // For now, using a simple orthographic projection (just discard z)
    return [point[0], point[1]];
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
    const { subvertexGap, dimension } = get(config);
    const subvertices = [];

    for (const [v1, v2] of edges) {
        const p1 = vertices[v1];
        const p2 = vertices[v2];

        // Calculate differences based on dimension
        const [dx, dy, dz] = dimension === 3
            ? [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]
            : [p2[0] - p1[0], p2[1] - p1[1], 0];

        const edgeLength = Math.sqrt(dx * dx + dy * dy + dz * dz);
        const subvertexCount = Math.max(1, Math.floor(edgeLength / subvertexGap));

        // Generate subvertices along the edge
        for (let i = 1; i <= subvertexCount; i++) {
            const t = i / (subvertexCount + 1);
            const position = dimension === 3
                ? [p1[0] + t * dx, p1[1] + t * dy, p1[2] + t * dz]
                : [p1[0] + t * dx, p1[1] + t * dy];

            subvertices.push({ position, edge: [v1, v2] });
        }
    }

    return subvertices;
}

export function generateRandomGraph(width, height) {
    const numVertices = Math.floor(Math.random() * 10) + 5;
    const vertices = [];
    const edges = [];
    const dimension = get(config).dimension;
    const defaultZ = get(config).defaultZ || 0;
    const zRange = 100; // Range for random Z values

    for (let i = 0; i < numVertices; i++) {
        if (dimension === 3) {
            const z = defaultZ + (Math.random() - 0.5) * zRange;
            vertices.push([Math.random() * width, Math.random() * height, z]);
        } else {
            vertices.push([Math.random() * width, Math.random() * height]);
        }
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
    const dimension = get(config).dimension;
    const defaultZ = get(config).defaultZ || 0;

    for (let i = 0; i < 3; i++) {
        if (dimension === 3) {
            const z = defaultZ + (i - 1) * 20; // Add some variance in Z
            vertices.push([width * 0.3, height * (0.2 + i * 0.3), z]);
        } else {
            vertices.push([width * 0.3, height * (0.2 + i * 0.3)]);
        }
    }

    for (let i = 0; i < 3; i++) {
        if (dimension === 3) {
            const z = defaultZ + (i - 1) * 20; // Add some variance in Z
            vertices.push([width * 0.7, height * (0.2 + i * 0.3), z]);
        } else {
            vertices.push([width * 0.7, height * (0.2 + i * 0.3)]);
        }
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
    const dimension = get(config).dimension;
    const defaultZ = get(config).defaultZ || 0;

    for (let i = 0; i < 2; i++) {
        if (dimension === 3) {
            const z = defaultZ + (i - 0.5) * 30; // Add some variance in Z
            vertices.push([width * 0.3, height * (0.3 + i * 0.4), z]);
        } else {
            vertices.push([width * 0.3, height * (0.3 + i * 0.4)]);
        }
    }

    for (let i = 0; i < 3; i++) {
        if (dimension === 3) {
            const z = defaultZ + (i - 1) * 20; // Add some variance in Z
            vertices.push([width * 0.7, height * (0.2 + i * 0.3), z]);
        } else {
            vertices.push([width * 0.7, height * (0.2 + i * 0.3)]);
        }
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
    const vertices = [];
    const dimension = get(config).dimension;
    const defaultZ = get(config).defaultZ || 0;
    
    if (dimension === 3) {
        vertices.push(
            [width * 0.3, height * 0.3, defaultZ - 20],
            [width * 0.7, height * 0.3, defaultZ + 20],
            [width * 0.7, height * 0.7, defaultZ + 20],
            [width * 0.3, height * 0.7, defaultZ - 20]
        );
    } else {
        vertices.push(
            [width * 0.3, height * 0.3],
            [width * 0.7, height * 0.3],
            [width * 0.7, height * 0.7],
            [width * 0.3, height * 0.7]
        );
    }

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
    const vertices = [];
    const dimension = get(config).dimension;
    const defaultZ = get(config).defaultZ || 0;
    
    if (dimension === 3) {
        vertices.push(
            [width * 0.3, height * 0.3, defaultZ - 30],
            [width * 0.7, height * 0.3, defaultZ + 30],
            [width * 0.7, height * 0.7, defaultZ + 30],
            [width * 0.3, height * 0.7, defaultZ - 30]
        );
    } else {
        vertices.push(
            [width * 0.3, height * 0.3],
            [width * 0.7, height * 0.3],
            [width * 0.7, height * 0.7],
            [width * 0.3, height * 0.7]
        );
    }

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