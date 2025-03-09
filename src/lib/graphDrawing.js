// src/lib/graphDrawing.js
import * as math from 'mathjs';
import { drawArrow, project3Dto2D } from '$lib/graphUtils';
import { calculateDifferential } from '$lib/energyCalculations';
import { canvasTransform, subvertices, config } from '$lib/stores';
import { get } from 'svelte/store';

export function drawGraph(
	ctx,
	width,
	height,
	vertices,
	edges,
	edgeProps,
	kernelMatrix,
	alpha,
	beta,
	disjointPairs
) {
	const { offsetX, offsetY, zoom } = get(canvasTransform);
	const subs = get(subvertices);
    const useSubverticesInEnergy = get(config).useSubverticesInEnergy;
    const dimension = get(config).dimension;

	ctx.save();
	ctx.clearRect(0, 0, width, height);

	ctx.scale(zoom, zoom);
	ctx.translate(offsetX / zoom, offsetY / zoom);

	drawEdges(ctx, vertices, edges, kernelMatrix, dimension);
	drawVertices(ctx, vertices, edges, alpha, beta, disjointPairs, dimension);
	drawMidpoints(ctx, edges, edgeProps, dimension);
	drawSubvertices(ctx, subs, useSubverticesInEnergy, dimension);
    
    // Draw energy inclusion status
    ctx.save();
    ctx.setTransform(1, 0, 0, 1, 0, 0); // Reset transform for text
    ctx.font = '14px Arial';
    ctx.fillStyle = 'black';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    ctx.fillText(
        `Mode: ${dimension}D | Energy includes ${useSubverticesInEnergy ? 'both vertices and subvertices' : 'only vertices'}`, 
        10, 10
    );
    ctx.restore();

	ctx.restore();
}

function drawEdges(ctx, vertices, edges, kernelMatrix, dimension) {
	const { zoom } = get(canvasTransform);

	if (kernelMatrix && math.isMatrix(kernelMatrix) && kernelMatrix.size()[0] === edges.length) {
		const maxKernelValue = math.max(kernelMatrix) || 1;

		edges.forEach((edge, i) => {
			const totalKernel = edges.reduce((sum, _, j) => {
				return i !== j ? sum + kernelMatrix.get([i, j]) : sum;
			}, 0);
			const avgKernel = totalKernel / (edges.length - 1 || 1);
			const normalizedValue = avgKernel / maxKernelValue;

			const blue = Math.round(255 * (1 - normalizedValue));
			const red = Math.round(255 * normalizedValue);

			ctx.strokeStyle = `rgb(${red}, 0, ${blue})`;
			ctx.lineWidth = (1 + normalizedValue * 4) / zoom;

			// Draw edges with proper projection if in 3D mode
			const v1 = vertices[edge[0]];
			const v2 = vertices[edge[1]];
			
			let p1, p2;
			if (dimension === 3) {
				p1 = project3Dto2D(v1);
				p2 = project3Dto2D(v2);
			} else {
				p1 = v1;
				p2 = v2;
			}

			ctx.beginPath();
			ctx.moveTo(p1[0], p1[1]);
			ctx.lineTo(p2[0], p2[1]);
			ctx.stroke();
			
			// If 3D, add depth cueing by adjusting opacity based on average Z
			if (dimension === 3) {
				const avgZ = (v1[2] + v2[2]) / 2;
				const zMin = -100; // Adjust these based on your scale
				const zMax = 100;
				// Make edges with higher Z more opaque (closer to viewer)
				const opacity = 0.4 + 0.6 * ((avgZ - zMin) / (zMax - zMin));
				ctx.globalAlpha = Math.max(0.1, Math.min(1, opacity));
			}
		});
		// Reset global alpha
		ctx.globalAlpha = 1.0;
	} else {
		ctx.strokeStyle = 'black';
		ctx.lineWidth = 1 / zoom;
		edges.forEach((edge) => {
			const v1 = vertices[edge[0]];
			const v2 = vertices[edge[1]];
			
			let p1, p2;
			if (dimension === 3) {
				p1 = project3Dto2D(v1);
				p2 = project3Dto2D(v2);
				
				// Add depth cueing
				const avgZ = (v1[2] + v2[2]) / 2;
				const zMin = -100;
				const zMax = 100;
				const opacity = 0.4 + 0.6 * ((avgZ - zMin) / (zMax - zMin));
				ctx.globalAlpha = Math.max(0.1, Math.min(1, opacity));
				
				// Optionally vary line width based on depth
				const depthScale = 0.5 + ((avgZ - zMin) / (zMax - zMin));
				ctx.lineWidth = (1 * depthScale) / zoom;
			} else {
				p1 = v1;
				p2 = v2;
			}
			
			ctx.beginPath();
			ctx.moveTo(p1[0], p1[1]);
			ctx.lineTo(p2[0], p2[1]);
			ctx.stroke();
		});
		// Reset global alpha
		ctx.globalAlpha = 1.0;
	}
}


function drawVertices(ctx, vertices, edges, alpha, beta, disjointPairs, dimension) {
    const { offsetX, offsetY, zoom } = get(canvasTransform);
    const gradient = calculateDifferential(vertices, edges, alpha, beta, disjointPairs);

    vertices.forEach((vertex, i) => {
        // Project from 3D to 2D if necessary
        const projectedPoint = dimension === 3 ? project3Dto2D(vertex) : vertex;
        
        // Calculate radius based on Z-value for 3D effect
        let radius = 5 / zoom;
        if (dimension === 3) {
            // Use absolute value mapping to ensure positive radius
            const zMin = -100;
            const zMax = 100;
            
            // Clamp Z value to our range
            const clampedZ = Math.max(zMin, Math.min(zMax, vertex[2]));
            
            // Calculate normalized position (0-1) within our range
            const zNormalized = (clampedZ - zMin) / (zMax - zMin);
            
            // Ensure radius is always positive with a minimum size
            radius = Math.max(3 / zoom, (3 + 2 * zNormalized) / zoom);
            
            // Adjust color based on depth
            const depthColor = Math.round(128 + 127 * zNormalized);
            ctx.fillStyle = `rgb(${depthColor}, ${depthColor}, 255)`;
        } else {
            ctx.fillStyle = 'blue';
        }
        
        ctx.beginPath();
        ctx.arc(projectedPoint[0], projectedPoint[1], radius, 0, 2 * Math.PI);
        ctx.fill();

        // Draw label at projected position
        const screenX = projectedPoint[0] * zoom + offsetX;
        const screenY = projectedPoint[1] * zoom + offsetY;

        ctx.save();
        ctx.setTransform(1, 0, 0, 1, 0, 0);

        ctx.fillStyle = 'black';
        ctx.font = '12px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        
        // In 3D mode, add Z coordinate to label
        if (dimension === 3) {
            ctx.fillText(`${i} (z:${vertex[2].toFixed(1)})`, screenX, screenY - 10);
        } else {
            ctx.fillText(i.toString(), screenX, screenY - 10);
        }

        // Draw gradient direction arrow
        const gradComponents = gradient[i];
        const gradX = -gradComponents[0];
        const gradY = -gradComponents[1];
        
        // For 3D, project gradient to 2D plane for visualization
        let dirX, dirY;
        if (dimension === 3) {
            const gradZ = -gradComponents[2];
            // Simple projection for now
            dirX = gradX;
            dirY = gradY;
            // Add indicator for Z component of gradient
            ctx.fillText(`∇z:${gradZ.toFixed(2)}`, screenX, screenY + 15);
        } else {
            dirX = gradX;
            dirY = gradY;
        }
        
        const magnitude = Math.sqrt(dirX * dirX + dirY * dirY) || 1e-6;
        const screenGradX = dirX * zoom;
        const screenGradY = dirY * zoom;
        const screenMagnitude = Math.sqrt(screenGradX * screenGradX + screenGradY * screenGradY);
        const arrowLength = 20;
        const arrowDirX = screenGradX / screenMagnitude;
        const arrowDirY = screenGradY / screenMagnitude;
        
        ctx.strokeStyle = 'purple';
        ctx.lineWidth = 2;
        drawArrow(ctx, screenX, screenY, arrowDirX, arrowDirY, arrowLength);

        ctx.restore();
    });
}


function drawMidpoints(ctx, edges, edgeProps, dimension) {
    const { offsetX, offsetY, zoom } = get(canvasTransform);

    if (!edgeProps || !edgeProps.edgeMidpoints || edgeProps.edgeMidpoints.length !== edges.length) {
        return;
    }

    edges.forEach((edge, i) => {
        const midpoint = edgeProps.edgeMidpoints[i];
        const length = edgeProps.edgeLengths[i];
        const tangent = edgeProps.edgeTangents[i];

        if (!midpoint || !length || !tangent) return;

        // Project midpoint if in 3D mode
        const projectedMidpoint = dimension === 3 ? project3Dto2D(midpoint) : midpoint;
        
        ctx.fillStyle = 'red';
        ctx.beginPath();
        ctx.arc(projectedMidpoint[0], projectedMidpoint[1], 2 / zoom, 0, 2 * Math.PI);
        ctx.fill();

        const screenX = projectedMidpoint[0] * zoom + offsetX;
        const screenY = projectedMidpoint[1] * zoom + offsetY;

        ctx.save();
        ctx.setTransform(1, 0, 0, 1, 0, 0);

        ctx.fillStyle = 'blue';
        ctx.font = '10px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        ctx.fillText(`${i}L ${length.toFixed(2)}`, screenX, screenY - 15);
        
        // Show tangent differently based on dimension
        if (dimension === 3) {
            ctx.fillText(`T [${tangent.map((t) => t.toFixed(1)).join(',')}]`, screenX, screenY - 5);
        } else {
            ctx.fillText(`T ${tangent.map((t) => t.toFixed(1))}`, screenX, screenY - 5);
        }
        
        ctx.fillText(`${edge[0]}, ${edge[1]}`, screenX, screenY + 15);

        // Draw tangent direction arrow
        ctx.strokeStyle = 'green';
        ctx.lineWidth = 1.5;
        
        // Project tangent to 2D if needed
        let projectedTangent;
        if (dimension === 3) {
            // Simple projection for now - just use X,Y components
            projectedTangent = [tangent[0], tangent[1]];
        } else {
            projectedTangent = tangent;
        }
        
        const tangentScreenX = projectedTangent[0] * zoom;
        const tangentScreenY = projectedTangent[1] * zoom;
        const tangentMagnitude = Math.sqrt(tangentScreenX * tangentScreenX + tangentScreenY * tangentScreenY) || 1;
        const dirX = tangentScreenX / tangentMagnitude;
        const dirY = tangentScreenY / tangentMagnitude;
        drawArrow(ctx, screenX, screenY, dirX, dirY, 20);

        ctx.restore();
    });
}


function drawSubvertices(ctx, subvertices, useInEnergy, dimension) {
    const { zoom } = get(canvasTransform);

    for (const subvertex of subvertices) {
        // Project subvertex position if in 3D
        const position = subvertex.position;
        const projectedPos = dimension === 3 ? project3Dto2D(position) : position;
        
        // Use different colors for subvertices based on whether they're included in energy
        ctx.fillStyle = useInEnergy ? 'purple' : 'black';
        ctx.strokeStyle = useInEnergy ? 'blue' : 'black';
        
        // Adjust size based on Z-coordinate if in 3D
        let radius = 4 / zoom;
        if (dimension === 3 && position.length >= 3) {
            const zMin = -100;
            const zMax = 100;
            const zNormalized = (position[2] - zMin) / (zMax - zMin);
            const sigmoid = (x) => 1 / (1 + Math.exp(-x / 4)); // TODO: softness of sigmoid based on the range of the z-vals
            const cappedZ = sigmoid(zNormalized * 6 - 3); // Scale and shift to center sigmoid around 0
            radius *= Math.max(0.5, cappedZ); // Ensure radius is always positive
        }
        
        ctx.beginPath();
        ctx.arc(projectedPos[0], projectedPos[1], radius, 0, 2 * Math.PI);
        ctx.fill();
        
        if (useInEnergy) {
            // Add a halo around subvertices to indicate they're included in energy
            ctx.lineWidth = 1 / zoom;
            ctx.beginPath();
            ctx.arc(projectedPos[0], projectedPos[1], radius + (3 / zoom), 0, 2 * Math.PI);
            ctx.stroke();
        }
    }
}


export function drawKernelMatrix(kernelCanvas, kernelMatrix) {
	if (!kernelCanvas || !kernelMatrix || !math.isMatrix(kernelMatrix)) return;

	const ctx = kernelCanvas.getContext('2d');
	const size = kernelMatrix.size();
	const padding = 50;
	const maxWidth = Math.min(window.innerWidth / 2, 800);
	const maxContentWidth = maxWidth - padding * 2;

	const cellSize = Math.min(
		50,
		Math.floor(maxContentWidth / size[1]),
		Math.floor(maxContentWidth / size[0])
	);

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