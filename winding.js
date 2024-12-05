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
