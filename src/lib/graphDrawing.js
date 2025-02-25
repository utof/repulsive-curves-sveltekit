import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls';
import * as math from 'mathjs';
import { calculateDifferential } from '$lib/energyCalculations';
import { drawArrow } from '$lib/graphUtils';
import { canvasTransform, subvertices, config } from '$lib/stores';
import { get } from 'svelte/store';

export function drawGraph(
    ctx2D,
    width,
    height,
    vertices,
    edges,
    edgeProps,
    kernelMatrix,
    alpha,
    beta,
    disjointPairs,
    threeRenderer,
    threeScene,
    threeCamera,
    threeControls
) {
    const dim = get(config).dim;
    if (dim === 2) {
        drawGraph2D(ctx2D, width, height, vertices, edges, edgeProps, kernelMatrix, alpha, beta, disjointPairs);
    } else if (dim === 3) {
        drawGraph3D(vertices, edges, threeRenderer, threeScene, threeCamera, threeControls);
    }
}

function drawGraph2D(ctx, width, height, vertices, edges, edgeProps, kernelMatrix, alpha, beta, disjointPairs) {
    const { offsetX, offsetY, zoom } = get(canvasTransform);
    const subs = get(subvertices);

    ctx.save();
    ctx.clearRect(0, 0, width, height);
    ctx.scale(zoom, zoom);
    ctx.translate(offsetX / zoom, offsetY / zoom);

    drawEdges2D(ctx, vertices, edges, kernelMatrix);
    drawVertices2D(ctx, vertices, edges, alpha, beta, disjointPairs);
    drawMidpoints2D(ctx, edges, edgeProps);
    drawSubvertices2D(ctx, subs);

    ctx.restore();
}

function drawEdges2D(ctx, vertices, edges, kernelMatrix) {
    const { zoom } = get(canvasTransform);

    if (kernelMatrix && math.isMatrix(kernelMatrix) && kernelMatrix.size()[0] === edges.length) {
        const maxKernelValue = math.max(kernelMatrix) || 1;

        edges.forEach((edge, i) => {
            const totalKernel = edges.reduce((sum, _, j) => i !== j ? sum + kernelMatrix.get([i, j]) : sum, 0);
            const avgKernel = totalKernel / (edges.length - 1 || 1);
            const normalizedValue = avgKernel / maxKernelValue;

            const blue = Math.round(255 * (1 - normalizedValue));
            const red = Math.round(255 * normalizedValue);

            ctx.strokeStyle = `rgb(${red}, 0, ${blue})`;
            ctx.lineWidth = (1 + normalizedValue * 4) / zoom;

            ctx.beginPath();
            ctx.moveTo(vertices[edge[0]][0], vertices[edge[0]][1]);
            ctx.lineTo(vertices[edge[1]][0], vertices[edge[1]][1]);
            ctx.stroke();
        });
    } else {
        ctx.strokeStyle = 'black';
        ctx.lineWidth = 1 / zoom;
        edges.forEach(edge => {
            ctx.beginPath();
            ctx.moveTo(vertices[edge[0]][0], vertices[edge[0]][1]);
            ctx.lineTo(vertices[edge[1]][0], vertices[edge[1]][1]);
            ctx.stroke();
        });
    }
}

function drawVertices2D(ctx, vertices, edges, alpha, beta, disjointPairs) {
    const { offsetX, offsetY, zoom } = get(canvasTransform);
    const gradient = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);

    vertices.forEach((vertex, i) => {
        ctx.beginPath();
        ctx.arc(vertex[0], vertex[1], 5 / zoom, 0, 2 * Math.PI);
        ctx.fillStyle = 'blue';
        ctx.fill();

        const screenX = vertex[0] * zoom + offsetX;
        const screenY = vertex[1] * zoom + offsetY;

        ctx.save();
        ctx.setTransform(1, 0, 0, 1, 0, 0);

        ctx.fillStyle = 'black';
        ctx.font = '12px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        ctx.fillText(i.toString(), screenX, screenY - 10);

        const gradX = -gradient[i][0];
        const gradY = -gradient[i][1];
        const magnitude = Math.sqrt(gradX * gradX + gradY * gradY) || 1e-6;
        const screenGradX = gradX * zoom;
        const screenGradY = gradY * zoom;
        const screenMagnitude = Math.sqrt(screenGradX * screenGradX + screenGradY * screenGradY);
        const arrowLength = 20;
        const dirX = screenGradX / screenMagnitude;
        const dirY = screenGradY / screenMagnitude;
        ctx.strokeStyle = 'purple';
        ctx.lineWidth = 2;
        drawArrow(ctx, screenX, screenY, dirX, dirY, arrowLength);

        ctx.restore();
    });
}

function drawMidpoints2D(ctx, edges, edgeProps) {
    const { offsetX, offsetY, zoom } = get(canvasTransform);

    if (!edgeProps || !edgeProps.edgeMidpoints || edgeProps.edgeMidpoints.length !== edges.length) return;

    edges.forEach((edge, i) => {
        const midpoint = edgeProps.edgeMidpoints[i];
        const length = edgeProps.edgeLengths[i];
        const tangent = edgeProps.edgeTangents[i];

        if (!midpoint || !length || !tangent) return;

        ctx.fillStyle = 'red';
        ctx.beginPath();
        ctx.arc(midpoint[0], midpoint[1], 2 / zoom, 0, 2 * Math.PI);
        ctx.fill();

        const screenX = midpoint[0] * zoom + offsetX;
        const screenY = midpoint[1] * zoom + offsetY;

        ctx.save();
        ctx.setTransform(1, 0, 0, 1, 0, 0);

        ctx.fillStyle = 'blue';
        ctx.font = '10px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        ctx.fillText(`${i}L ${length.toFixed(6)}`, screenX, screenY - 15);
        ctx.fillText(`T ${tangent.map(t => t.toFixed(1))}`, screenX, screenY - 5);
        ctx.fillText(`${edge[0]}, ${edge[1]}`, screenX, screenY + 15);

        ctx.strokeStyle = 'green';
        ctx.lineWidth = 1.5;
        const tangentScreenX = tangent[0] * zoom;
        const tangentScreenY = tangent[1] * zoom;
        const tangentMagnitude = Math.sqrt(tangentScreenX * tangentScreenX + tangentScreenY * tangentScreenY) || 1;
        const dirX = tangentScreenX / tangentMagnitude;
        const dirY = tangentScreenY / tangentMagnitude;
        drawArrow(ctx, screenX, screenY, dirX, dirY, 20);

        ctx.restore();
    });
}

function drawSubvertices2D(ctx, subvertices) {
    const { zoom } = get(canvasTransform);

    ctx.fillStyle = 'black';
    for (const subvertex of subvertices) {
        const [x, y] = subvertex.position;
        ctx.beginPath();
        ctx.arc(x, y, 4 / zoom, 0, 2 * Math.PI);
        ctx.fill();
    }
}

function drawGraph3D(vertices, edges, renderer, scene, camera, controls) {
    while (scene.children.length > 0) {
        scene.remove(scene.children[0]);
    }

    const ambientLight = new THREE.AmbientLight(0xffffff, 0.5);
    scene.add(ambientLight);

    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.5);
    directionalLight.position.set(0, 1, 1).normalize();
    scene.add(directionalLight);

    vertices.forEach(vertex => {
        const geometry = new THREE.SphereGeometry(5, 16, 16);
        const material = new THREE.MeshPhongMaterial({ color: 0xffffff });
        const sphere = new THREE.Mesh(geometry, material);
        sphere.position.set(vertex[0], vertex[1], vertex[2]);
        scene.add(sphere);
    });

    edges.forEach(([v1, v2]) => {
        const material = new THREE.LineBasicMaterial({ color: 0xff0000, linewidth: 5 }); // Thicker lines
        const points = [
            new THREE.Vector3(vertices[v1][0], vertices[v1][1], vertices[v1][2]),
            new THREE.Vector3(vertices[v2][0], vertices[v2][1], vertices[v2][2])
        ];
        const geometry = new THREE.BufferGeometry().setFromPoints(points);
        const line = new THREE.Line(geometry, material);
        scene.add(line);
    });

    renderer.render(scene, camera);
    controls.update();
}

export function drawKernelMatrix(kernelCanvas, kernelMatrix) {
    // Existing implementation remains unchanged
    if (!kernelCanvas || !kernelMatrix || !math.isMatrix(kernelMatrix)) return;

    const ctx = kernelCanvas.getContext('2d');
    const size = kernelMatrix.size();
    const padding = 50;
    const maxWidth = Math.min(window.innerWidth / 2, 800);
    const maxContentWidth = maxWidth - padding * 2;

    const cellSize = Math.min(50, Math.floor(maxContentWidth / size[1]), Math.floor(maxContentWidth / size[0]));

    const canvasWidth = Math.min(maxWidth, padding * 2 + size[1] * cellSize);
    const canvasHeight = padding * 2 + size[0] * cellSize;

    kernelCanvas.width = canvasWidth;
    kernelCanvas.height = canvasHeight;
    kernelCanvas.style.width = `${canvasWidth}px`;
    kernelCanvas.style.height = `${canvasHeight}px`;

    ctx.clearRect(0, 0, canvasWidth, canvasHeight);

    const maxVal = math.max(kernelMatrix) || 1;

    for (let i = 0; i < size[0]; i++) {
        for (let j = 0; j < size[1]; j++) {
            const value = kernelMatrix.get([i, j]);
            const intensity = value / maxVal;
            const r = Math.round(255 - intensity * 255);
            const g = Math.round(255 - intensity * 255);
            const b = 255;
            const minOpacity = 0.1;
            const opacity = minOpacity + (1 - minOpacity) * intensity;
            ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${opacity})`;
            ctx.fillRect(padding + j * cellSize, padding + i * cellSize, cellSize, cellSize);
        }
    }

    ctx.strokeStyle = 'black';
    ctx.lineWidth = 1;
    for (let i = 0; i <= size[0]; i++) {
        ctx.beginPath();
        ctx.moveTo(padding, padding + i * cellSize);
        ctx.lineTo(padding + size[1] * cellSize, padding + i * cellSize);
        ctx.stroke();
    }
    for (let j = 0; j <= size[1]; j++) {
        ctx.beginPath();
        ctx.moveTo(padding + j * cellSize, padding);
        ctx.lineTo(padding + j * cellSize, padding + size[0] * cellSize);
        ctx.stroke();
    }

    ctx.fillStyle = 'black';
    ctx.font = '12px Arial';
    ctx.textAlign = 'center';
    for (let j = 0; j < size[1]; j++) {
        ctx.fillText(j.toString(), padding + j * cellSize + cellSize / 2, padding - 10);
    }

    ctx.textAlign = 'right';
    for (let i = 0; i < size[0]; i++) {
        ctx.fillText(i.toString(), padding - 5, padding + i * cellSize + cellSize / 2);
    }
}