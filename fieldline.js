// Drawing of field lines and contours.
//
// Field lines are only ever drawn for the analytic 2D models (see fieldmodels.js and
// https://github.com/apontzen/maxwell/issues/3): their geometry comes from fieldlines2d.js
// and is drawn here, together with direction arrows and Gaussian-surface crossings.
// Potential contours of the 3D point-charge model (explorer, quiver mode) still use the
// contour tracer in the WASM module.

import { generate_potential_contours_at_levels } from './maxwell/out/maxwell.js';
import { drawGaussianSurfaceCrossings } from './gaussian_surface.js';
import { electricFieldLines, magneticFieldLines, potentialContours } from './fieldlines2d.js';

export function drawDirectionArrow(x_position, y_position, u, v, ctx, color) {
    ctx.save();
    if (color !== undefined) ctx.strokeStyle = color;
    ctx.lineWidth = 1;
    ctx.translate(x_position, y_position);
    ctx.rotate(Math.atan2(v, u));
    ctx.beginPath();
    ctx.moveTo(-8, -5);
    ctx.lineTo(0, 0);
    ctx.lineTo(-8, 5);
    ctx.stroke();
    ctx.restore();
}

function strokePolyline(ctx, points, closed = false) {
    if (points.length < 2) return;
    ctx.beginPath();
    ctx.moveTo(points[0][0], points[0][1]);
    for (let i = 1; i < points.length; i++) ctx.lineTo(points[i][0], points[i][1]);
    if (closed) ctx.closePath();
    ctx.stroke();
}

/** Maximal runs of consecutive points inside the canvas. */
function visibleRuns(points, rect) {
    const runs = [];
    let run = null;
    for (const p of points) {
        const inside = p[0] >= 0 && p[1] >= 0 && p[0] <= rect.width && p[1] <= rect.height;
        if (inside) {
            if (!run) { run = []; runs.push(run); }
            run.push(p);
        } else {
            run = null;
        }
    }
    return runs.filter(r => r.length > 1);
}

/** The point half-way along a polyline, measured by arc length. */
function arcLengthMidpoint(points) {
    let total = 0;
    for (let i = 1; i < points.length; i++)
        total += Math.hypot(points[i][0] - points[i - 1][0], points[i][1] - points[i - 1][1]);
    let target = total / 2;
    for (let i = 1; i < points.length; i++) {
        const seg = Math.hypot(points[i][0] - points[i - 1][0], points[i][1] - points[i - 1][1]);
        if (seg >= target && seg > 0) {
            const t = target / seg;
            return [points[i - 1][0] + t * (points[i][0] - points[i - 1][0]),
                    points[i - 1][1] + t * (points[i][1] - points[i - 1][1])];
        }
        target -= seg;
    }
    return points[points.length - 1];
}

function pathLength(points) {
    let total = 0;
    for (let i = 1; i < points.length; i++)
        total += Math.hypot(points[i][0] - points[i - 1][0], points[i][1] - points[i - 1][1]);
    return total;
}

/**
 * Candidate direction-arrow positions for a polyline: the middle of each visible stretch. A
 * closed loop that is entirely visible gets two, at its leftmost and rightmost points, which
 * move continuously as charges are dragged. Each candidate is weighted by the length of the
 * stretch it belongs to.
 */
function arrowCandidates(points, closed, rect) {
    const runs = visibleRuns(points, rect);
    if (closed && runs.length === 1 && runs[0].length === points.length) {
        let left = points[0], right = points[0];
        for (const p of points) {
            if (p[0] < left[0]) left = p;
            if (p[0] > right[0]) right = p;
        }
        const weight = pathLength(points);
        return [{ p: left, weight }, { p: right, weight }];
    }
    return runs.map(run => ({ p: arcLengthMidpoint(run), weight: pathLength(run) }));
}

const MIN_ARROW_SEPARATION = 18; // px
const MIN_ARROW_DISTANCE_FROM_SOURCE = 16; // px (the drawn charges have radius 10)

/**
 * Draw arrows at the candidate positions, longest lines first, skipping any that would crowd an
 * arrow already placed or sit on top of a charge (e.g. on the many small loops around a wire).
 */
function drawArrows(ctx, model, candidates, color) {
    candidates.sort((a, b) => b.weight - a.weight);
    const placed = [];
    for (const { p } of candidates) {
        if (model.sources.some(s => Math.hypot(p[0] - s.x, p[1] - s.y) < MIN_ARROW_DISTANCE_FROM_SOURCE)) continue;
        if (placed.some(q => Math.hypot(p[0] - q[0], p[1] - q[1]) < MIN_ARROW_SEPARATION)) continue;
        placed.push(p);
        const { u, v } = model.fieldAt(p[0], p[1]);
        drawDirectionArrow(p[0], p[1], u, v, ctx, color);
    }
}

function crossings(ctx, model, points, closed, gaussianSurfacePoints) {
    if (gaussianSurfacePoints.length === 0) return;
    const pts = points.map(([x, y]) => ({ x, y }));
    if (closed && pts.length > 1) pts.push(pts[0]);
    drawGaussianSurfaceCrossings(pts, gaussianSurfacePoints, model.fieldToBuffer, null, ctx);
}

/** Electric field lines of a 2D line-charge model. */
export function drawElectricFieldLines2D(ctx, model, rect, gaussianSurfacePoints) {
    const color = gaussianSurfacePoints.length > 0 ? 'grey' : 'black';
    const lines = electricFieldLines(model, { width: rect.width, height: rect.height });
    const candidates = [];
    for (const line of lines) {
        ctx.save();
        ctx.strokeStyle = color;
        ctx.lineWidth = 1.0;
        strokePolyline(ctx, line.points);
        ctx.restore();
        candidates.push(...arrowCandidates(line.points, false, rect));
    }
    drawArrows(ctx, model, candidates, color);
    for (const line of lines) crossings(ctx, model, line.points, false, gaussianSurfacePoints);
}

/** Magnetic field lines (contours of A_z) of a 2D straight-wire model. */
export function drawMagneticFieldLines2D(ctx, model, rect, gaussianSurfacePoints) {
    const lines = magneticFieldLines(model, { width: rect.width, height: rect.height });
    const candidates = [];
    for (const line of lines) {
        ctx.save();
        ctx.strokeStyle = 'black';
        ctx.lineWidth = 1.0;
        strokePolyline(ctx, line.points, line.closed);
        ctx.restore();
        candidates.push(...arrowCandidates(line.points, line.closed, rect));
    }
    drawArrows(ctx, model, candidates, 'black');
    for (const line of lines) crossings(ctx, model, line.points, line.closed, gaussianSurfacePoints);
}

/** Equipotentials of a 2D line-charge model: positive red, negative blue, zero grey. */
export function drawPotentialContours2D(ctx, model, rect) {
    const contours = potentialContours(model, { width: rect.width, height: rect.height });
    for (const c of contours) {
        ctx.save();
        ctx.strokeStyle = c.level > 0 ? 'red' : (c.level < 0 ? 'blue' : 'grey');
        ctx.lineWidth = 1.0;
        strokePolyline(ctx, c.points, c.closed);
        ctx.restore();
    }
}

/** Equipotentials of the 3D point-charge model, traced by the WASM module. */
export function drawPotentialContours(field, levels, ctx, color) {
    const contours = generate_potential_contours_at_levels(field, levels);
    for (let contour of contours) {
        ctx.strokeStyle = color;
        ctx.lineWidth = 1.0;
        ctx.beginPath();
        for (let i = 0; i < contour.length; i++) {
            const [x, y] = contour[i];
            if (i === 0) {
                ctx.moveTo(x, y);
            } else {
                ctx.lineTo(x, y);
            }
        }
        ctx.stroke();
    }
}
