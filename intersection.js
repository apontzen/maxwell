/**
 * Module to find the intersection points between two curves efficiently.
 *
 * @param {Array} c1 - First curve, array of points [{x, y}, {x, y}, ...]
 * @param {Array} c2 - Second curve, array of points [{x, y}, {x, y}, ...]
 * @returns {Array} - Array of intersection points [{x, y}, ...]
 */
export function getAllIntersections(c1, c2) {
    // Convert points to line segments with bounding boxes
    const segments1 = makeLineSegments(c1);
    const segments2 = makeLineSegments(c2);
  
    // Sort segments by minX
    segments1.sort((a, b) => a.bbox.minX - b.bbox.minX);
    segments2.sort((a, b) => a.bbox.minX - b.bbox.minX);
  
    const intersections = [];
  
    // For each segment in segments1, find potential overlapping segments in segments2
    let index2 = 0;
    const len2 = segments2.length;
  
    for (let i = 0; i < segments1.length; i++) {
      const s1 = segments1[i];
      const possibleSegments = [];
  
      // Move index2 to the first segment in segments2 that could overlap with s1
      while (
        index2 < len2 &&
        segments2[index2].bbox.maxX < s1.bbox.minX
      ) {
        index2++;
      }
  
      // Collect all segments in segments2 that could overlap with s1 in the x-axis
      let j = index2;
      while (j < len2 && segments2[j].bbox.minX <= s1.bbox.maxX) {
        const s2 = segments2[j];
  
        // Check for overlap in y-axis
        if (
          s1.bbox.maxY >= s2.bbox.minY &&
          s1.bbox.minY <= s2.bbox.maxY
        ) {
          // Check for intersection
          const intersection = getLineIntersection(s1, s2);
          if (intersection) {
            intersections.push(intersection);
          }
        }
        j++;
      }
    }
  
    return intersections;
  }
  
  // Helper function to create line segments from points with bounding boxes
  function makeLineSegments(points) {
    const segments = [];
    for (let i = 0; i < points.length - 1; i++) {
      const p1 = points[i];
      const p2 = points[i + 1];
      const minX = Math.min(p1.x, p2.x);
      const minY = Math.min(p1.y, p2.y);
      const maxX = Math.max(p1.x, p2.x);
      const maxY = Math.max(p1.y, p2.y);
      segments.push({
        x1: p1.x,
        y1: p1.y,
        x2: p2.x,
        y2: p2.y,
        bbox: { minX, minY, maxX, maxY },
      });
    }
    return segments;
  }
  
  // Helper function to check if two line segments intersect
  function getLineIntersection(l1, l2) {
    const { x1, y1, x2, y2 } = l1;
    const { x1: x3, y1: y3, x2: x4, y2: y4 } = l2;
  
    const denom =
      (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
  
    if (denom === 0) return null; // Parallel or colinear
  
    const ua =
      ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denom;
    const ub =
      ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denom;
  
    if (ua >= 0 && ua <= 1 && ub >= 0 && ub <= 1) {
      return {
        x: x1 + ua * (x2 - x1),
        y: y1 + ua * (y2 - y1),
      };
    }
    return null;
  }
