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
        return [{ p: left, weight, run: points, closed: true }, { p: right, weight, run: points, closed: true }];
    }
    return runs.map(run => ({ p: arcLengthMidpoint(run), weight: pathLength(run), run, closed: false }));
}

/** The point of a polyline closest to q, its distance, and its arc length from either end. */
function closestOnPolyline(run, q, closed) {
    let best = null, bestD = Infinity, s = 0, bestS = 0;
    const n = closed ? run.length : run.length - 1;
    for (let i = 0; i < n; i++) {
        const a = run[i], b = run[(i + 1) % run.length];
        const dx = b[0] - a[0], dy = b[1] - a[1], l2 = dx * dx + dy * dy, l = Math.sqrt(l2);
        const t = l2 > 0 ? Math.max(0, Math.min(1, ((q[0] - a[0]) * dx + (q[1] - a[1]) * dy) / l2)) : 0;
        const x = a[0] + t * dx, y = a[1] + t * dy, d = Math.hypot(q[0] - x, q[1] - y);
        if (d < bestD) { bestD = d; best = [x, y]; bestS = s + t * l; }
        s += l;
    }
    return { p: best, d: bestD, fromStart: bestS, fromEnd: s - bestS };
}

// Arrows that stay where they were (see drawArrows), per canvas
const previousArrows = new WeakMap();
const ARROW_STICK_DISTANCE = 3;   // px: a line that moved less than this keeps its arrow in place
const ARROW_STICK_END_MARGIN = 12; // px: ... unless the arrow would end up this close to an end of its stretch

const MIN_ARROW_SEPARATION = 18; // px
const MIN_ARROW_DISTANCE_FROM_SOURCE = 16; // px (the drawn charges have radius 10)

/**
 * Draw arrows at the candidate positions, longest lines first, skipping any that would crowd an
 * arrow already placed or sit on top of a charge (e.g. on the many small loops around a wire).
 *
 * With `memoryKey` (the dynamic solver, whose lines evolve from frame to frame), an arrow stays where
 * it was in the previous frame, moved onto its line, as long as the line has moved less than
 * ARROW_STICK_DISTANCE: otherwise the midpoint of a line that is unchanged except near its far end
 * (e.g. where it leaves the canvas) slides, and arrows on still lines would jitter.
 */
function drawArrows(ctx, model, candidates, color, memoryKey = null) {
    candidates.sort((a, b) => b.weight - a.weight);
    const previous = (memoryKey && previousArrows.get(memoryKey)) || [];
    const placed = [];
    for (const c of candidates) {
        let p = c.p;
        for (const q of previous) {
            const near = closestOnPolyline(c.run, q, c.closed);
            if (near.d < ARROW_STICK_DISTANCE && (c.closed ||
                (near.fromStart > ARROW_STICK_END_MARGIN && near.fromEnd > ARROW_STICK_END_MARGIN))) { p = near.p; break; }
        }
        if (model.sources.some(s => Math.hypot(p[0] - s.x, p[1] - s.y) < MIN_ARROW_DISTANCE_FROM_SOURCE)) continue;
        if (placed.some(q => Math.hypot(p[0] - q[0], p[1] - q[1]) < MIN_ARROW_SEPARATION)) continue;
        placed.push(p);
        const { u, v } = model.fieldAt(p[0], p[1]);
        drawDirectionArrow(p[0], p[1], u, v, ctx, color);
    }
    if (memoryKey) previousArrows.set(memoryKey, placed);
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
        strokePolyline(ctx, line.points, !!line.closed);
        ctx.restore();
        candidates.push(...arrowCandidates(line.points, !!line.closed, rect));
    }
    drawArrows(ctx, model, candidates, color, model.stickyArrows ? ctx.canvas : null);
    for (const line of lines) crossings(ctx, model, line.points, !!line.closed, gaussianSurfacePoints);
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
