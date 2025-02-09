import { getAllIntersections } from './intersection.js';
import { drawArrow } from './draw.js';

export function windingNumber(points, testPoint) {
    const { x: px, y: py } = testPoint;
    let windingNum = 0;

    for (let i = 0; i < points.length; i++) {
        // Get the current point and the next point, wrapping around for a closed curve
        const { x: x1, y: y1 } = points[i];
        const { x: x2, y: y2 } = points[(i + 1) % points.length];

        // Check if the segment crosses the horizontal line from the test point
        if (y1 <= py) {
            if (y2 > py) {
                // Check if the segment is counterclockwise
                const isLeft = (x2 - x1) * (py - y1) - (px - x1) * (y2 - y1);
                if (isLeft > 0) {
                    windingNum++;
                }
            }
        } else {
            if (y2 <= py) {
                // Check if the segment is clockwise
                const isLeft = (x2 - x1) * (py - y1) - (px - x1) * (y2 - y1);
                if (isLeft < 0) {
                    windingNum--;
                }
            }
        }
    }

    return windingNum;
}

export function drawGaussianSurfaceCrossings(fieldlinePoints, gaussianSurfacePoints, compute_field_to_buffer, field, ctx) {
    let intersections = getAllIntersections(fieldlinePoints, gaussianSurfacePoints);
    let buffer = new Float64Array(2);

    for (let intersection of intersections) {
        compute_field_to_buffer(field, intersection.x, intersection.y, buffer);
        let u = buffer[0];
        let v = buffer[1];
        const norm = Math.sqrt(u * u + v * v);
        u*=80.0/norm;
        v*=80.0/norm;

        const testHead = {x: intersection.x + u*0.01, y: intersection.y + v*0.01};
        const testTail = {x: intersection.x - u*0.01, y: intersection.y - v*0.01};

        let headWinding = windingNumber(gaussianSurfacePoints, testHead);
        let tailWinding = windingNumber(gaussianSurfacePoints, testTail);

        if (Math.abs(headWinding)>1 || Math.abs(tailWinding)>1)
            continue;

        let color;
        if(headWinding == 0)
            color = 'green';
        else   
            color = 'purple';
        drawArrow(ctx, intersection.x, intersection.y, u, v, color, 2);

    }
}