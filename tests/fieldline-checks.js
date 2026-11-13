// Pure, environment-independent checker functions for 2D field-line geometry.
//
// Nothing here touches the DOM, Playwright or the field-line implementation, so the
// module can be imported in Node (see fieldline-checks.spec.js, which validates every
// checker against exact analytic data) and in the browser (fieldline-rules.spec.js).
//
// Conventions
//   * A polyline is an array of [x, y] points. A "line" (as returned by
//     electricFieldLines) is {points, from, to, level}.
//   * A source is {x, y, charge}.
//   * Angles are radians unless the name says otherwise; angles of directions are
//     measured with atan2(dy, dx) in canvas coordinates (y points down; nothing here
//     depends on the handedness, it just has to be used consistently).
//   * A line "end" is 'start' (points[0]) or 'end' (points[points.length - 1]).

export const TAU = 2 * Math.PI;
export const DEG = Math.PI / 180;

// ---------------------------------------------------------------------------
// Small geometry helpers
// ---------------------------------------------------------------------------

/** Angle reduced to [0, 2 pi) */
export function normAngle(a) {
    const r = a % TAU;
    return r < 0 ? r + TAU : r;
}

/** Signed difference a - b reduced to (-pi, pi] */
export function angleDiff(a, b) {
    let d = normAngle(a - b);
    if (d > Math.PI) d -= TAU;
    return d;
}

export const dist = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1]);

export function endPoint(points, end) {
    return end === 'start' ? points[0] : points[points.length - 1];
}

/** Distance from point p to segment a-b */
export function pointSegmentDistance(p, a, b) {
    const dx = b[0] - a[0], dy = b[1] - a[1];
    const len2 = dx * dx + dy * dy;
    let t = len2 === 0 ? 0 : ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2;
    t = Math.max(0, Math.min(1, t));
    return Math.hypot(p[0] - (a[0] + t * dx), p[1] - (a[1] + t * dy));
}

/** Distance from p to a polyline (open) */
export function pointPolylineDistance(p, points) {
    if (points.length === 1) return dist(p, points[0]);
    let best = Infinity;
    for (let i = 0; i + 1 < points.length; i++)
        best = Math.min(best, pointSegmentDistance(p, points[i], points[i + 1]));
    return best;
}

/** Distance from p to the boundary of a closed polygon */
export function pointPolygonDistance(p, polygon) {
    let best = Infinity;
    for (let i = 0; i < polygon.length; i++)
        best = Math.min(best, pointSegmentDistance(p, polygon[i], polygon[(i + 1) % polygon.length]));
    return best;
}

/** Even-odd point-in-polygon test */
export function pointInPolygon(p, polygon) {
    let inside = false;
    for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
        const [xi, yi] = polygon[i], [xj, yj] = polygon[j];
        if ((yi > p[1]) !== (yj > p[1]) && p[0] < ((xj - xi) * (p[1] - yi)) / (yj - yi) + xi)
            inside = !inside;
    }
    return inside;
}

/** Total length of a polyline */
export function polylineLength(points) {
    let s = 0;
    for (let i = 0; i + 1 < points.length; i++) s += dist(points[i], points[i + 1]);
    return s;
}

// ---------------------------------------------------------------------------
// Attachment of lines to sources
// ---------------------------------------------------------------------------

/**
 * For each source, the line ends attached to it (within `tol` px of its position).
 * Returns an array (one entry per source) of [{line, end}] with `line` an index into `lines`.
 */
export function attachedEnds(lines, sources, tol = 1) {
    const out = sources.map(() => []);
    lines.forEach((line, li) => {
        for (const end of ['start', 'end']) {
            const p = endPoint(line.points, end);
            sources.forEach((s, si) => {
                if (Math.hypot(p[0] - s.x, p[1] - s.y) <= tol) out[si].push({ line: li, end });
            });
        }
    });
    return out;
}

/** Number of line ends attached to each source (an end counts if within `tol` px of the source). */
export function countAttachedEnds(lines, sources, tol = 1) {
    return attachedEnds(lines, sources, tol).map(a => a.length);
}

/**
 * Point on the polyline at Euclidean distance `d` from `centre`, found by walking away from
 * the given end of the polyline (first crossing of the circle of radius d). If the polyline never
 * gets that far, its far end is returned.
 */
export function pointAtDistance(points, end, centre, d) {
    const pts = end === 'start' ? points : [...points].reverse();
    const r = p => Math.hypot(p[0] - centre[0], p[1] - centre[1]);
    for (let i = 1; i < pts.length; i++) {
        const r0 = r(pts[i - 1]), r1 = r(pts[i]);
        if (r1 >= d && r0 < d) {
            const t = (d - r0) / (r1 - r0);
            return [pts[i - 1][0] + t * (pts[i][0] - pts[i - 1][0]), pts[i - 1][1] + t * (pts[i][1] - pts[i - 1][1])];
        }
    }
    return pts[pts.length - 1];
}

/**
 * Attachment angle of a line at a source: direction (from the source) of the point where the
 * polyline, followed away from the given end, first reaches distance `d` (default 15 px)
 * from the source. For d -> 0 this tends to the tangent direction of the line at the source;
 * for a circular arc of radius R it differs from it by asin(d / 2R). Returned in [0, 2 pi).
 */
export function attachmentAngle(points, end, source, d = 15) {
    const p = pointAtDistance(points, end, [source.x, source.y], d);
    return normAngle(Math.atan2(p[1] - source.y, p[0] - source.x));
}

/**
 * Sorted attachment angles of all line ends attached to each source.
 * Returns an array (per source) of ascending angles in [0, 2 pi).
 */
export function attachmentAngles(lines, sources, { d = 15, tol = 1 } = {}) {
    return attachedEnds(lines, sources, tol).map((ends, si) =>
        ends.map(({ line, end }) => attachmentAngle(lines[line].points, end, sources[si], d))
            .sort((a, b) => a - b));
}

/**
 * Compare two sets of angles that are only defined up to cyclic order (e.g. the attachment
 * angles at one charge in two frames). `offset` is added to every angle of A first (use it for a
 * rotation). Both sets are sorted and the best cyclic alignment is found; the result is
 * {ok: sizes equal, maxDiff: largest remaining difference in radians (Infinity if sizes differ),
 *  shift}. Differences are measured on the circle.
 */
export function matchAngleSets(A, B, offset = 0) {
    if (A.length !== B.length) return { ok: false, maxDiff: Infinity, shift: null, sizeA: A.length, sizeB: B.length };
    const n = A.length;
    if (n === 0) return { ok: true, maxDiff: 0, shift: 0 };
    const a = A.map(x => normAngle(x + offset)).sort((x, y) => x - y);
    const b = B.map(normAngle).sort((x, y) => x - y);
    let best = { maxDiff: Infinity, shift: 0 };
    for (let s = 0; s < n; s++) {
        let m = 0;
        for (let i = 0; i < n; i++) m = Math.max(m, Math.abs(angleDiff(a[i], b[(i + s) % n])));
        if (m < best.maxDiff) best = { maxDiff: m, shift: s };
    }
    return { ok: true, ...best };
}

// ---------------------------------------------------------------------------
// Endpoint consistency
// ---------------------------------------------------------------------------

/**
 * Check that every line end is either on a source (within `tol` px) or marked 'infinity' and
 * lies more than `farFactor` x the view diagonal away from the view centre, and that the
 * from / to labels agree with the geometry. Returns a list of problem strings (empty = fine).
 * `lines` must be {points, from, to}; for lines without from/to (magnetic) pass labelled=false
 * and then open ends must be far away and closed lines must be marked closed.
 */
export function checkEndpoints(lines, sources, view, { tol = 1, farFactor = 3, labelled = true, diskRadius = null } = {}) {
    const problems = [];
    const cx = view.width / 2, cy = view.height / 2;
    const far = farFactor * Math.hypot(view.width, view.height);
    lines.forEach((line, li) => {
        const name = `line ${li} (level ${line.level})`;
        if (!line.points || line.points.length < 2) {
            problems.push(`${name} has fewer than 2 points`);
            return;
        }
        for (const end of ['start', 'end']) {
            const label = end === 'start' ? line.from : line.to;
            const p = endPoint(line.points, end);
            let nearest = -1, nearestDist = Infinity;
            sources.forEach((s, si) => {
                const dd = Math.hypot(p[0] - s.x, p[1] - s.y);
                if (dd < nearestDist) { nearestDist = dd; nearest = si; }
            });
            const onSource = nearestDist <= tol;
            const farAway = Math.hypot(p[0] - cx, p[1] - cy) > far;
            const where = `${name} ${end} (${p[0].toFixed(1)}, ${p[1].toFixed(1)})`;
            if (!labelled) {
                if (line.closed) continue;
                // contours may be clipped where they enter the drawn disk of a source
                if (diskRadius !== null && Math.abs(nearestDist - diskRadius) <= tol) continue;
                if (!farAway) problems.push(`${where} is an open end but is not more than ${farFactor} view diagonals from the view centre`);
                continue;
            }
            if (label === 'infinity') {
                if (onSource) problems.push(`${where} is labelled 'infinity' but lies on source ${nearest}`);
                else if (!farAway) problems.push(`${where} is labelled 'infinity' but is only ${Math.hypot(p[0] - cx, p[1] - cy).toFixed(0)} px from the view centre (needs > ${far.toFixed(0)})`);
            } else if (typeof label === 'number') {
                if (label < 0 || label >= sources.length) problems.push(`${where} has out-of-range ${end === 'start' ? 'from' : 'to'} = ${label}`);
                else if (Math.hypot(p[0] - sources[label].x, p[1] - sources[label].y) > tol)
                    problems.push(`${where} is labelled ${end === 'start' ? 'from' : 'to'} = ${label} but is ${Math.hypot(p[0] - sources[label].x, p[1] - sources[label].y).toFixed(2)} px from that source`);
            } else {
                problems.push(`${where} has invalid ${end === 'start' ? 'from' : 'to'} label ${JSON.stringify(label)}`);
            }
        }
    });
    return problems;
}

// ---------------------------------------------------------------------------
// Crossings with a closed polygon (Gauss's law)
// ---------------------------------------------------------------------------

/** Twice the signed area (positive for counter-clockwise in the usual x right / y up sense) */
function signedArea2(polygon) {
    let a = 0;
    for (let i = 0; i < polygon.length; i++) {
        const [x0, y0] = polygon[i], [x1, y1] = polygon[(i + 1) % polygon.length];
        a += x0 * y1 - x1 * y0;
    }
    return a;
}

/**
 * Signed crossings of a polyline with a closed polygon. A crossing is +1 if the polyline leaves
 * the polygon (crosses outwards) and -1 if it enters. For a field line oriented from + to - the
 * sum over lines is the flux in units of 1/LINES_PER_UNIT_CHARGE, i.e. 4 x enclosed charge.
 * Every segment-edge intersection is counted (a segment may cross twice), so sparse polylines
 * are handled exactly. Returns {out, in, net}.
 */
export function signedCrossings(points, polygon) {
    const orient = signedArea2(polygon) > 0 ? 1 : -1;
    let out = 0, inn = 0;
    for (let i = 0; i + 1 < points.length; i++) {
        const p = points[i], q = points[i + 1];
        const rx = q[0] - p[0], ry = q[1] - p[1];
        for (let j = 0; j < polygon.length; j++) {
            const a = polygon[j], b = polygon[(j + 1) % polygon.length];
            const sx = b[0] - a[0], sy = b[1] - a[1];
            const denom = rx * sy - ry * sx;
            if (denom === 0) continue;
            const apx = a[0] - p[0], apy = a[1] - p[1];
            const t = (apx * sy - apy * sx) / denom;
            const u = (apx * ry - apy * rx) / denom;
            // half-open intervals so that a crossing through a vertex is counted exactly once
            if (t < 0 || t >= 1 || u < 0 || u >= 1) continue;
            // outward normal of the edge is (sy, -sx) for counter-clockwise polygons
            const outward = orient * (rx * sy - ry * sx) > 0;
            if (outward) out++; else inn++;
        }
    }
    return { out, in: inn, net: out - inn };
}

/** Net signed crossings of a set of lines (objects with .points) with a polygon. */
export function netCrossings(lines, polygon) {
    let net = 0, out = 0, inn = 0;
    for (const line of lines) {
        const c = signedCrossings(line.points, polygon);
        net += c.net; out += c.out; inn += c.in;
    }
    return { net, out, in: inn };
}

/** Sum of charges of the sources inside a polygon */
export function enclosedCharge(sources, polygon) {
    return sources.filter(s => pointInPolygon([s.x, s.y], polygon)).reduce((t, s) => t + s.charge, 0);
}

// ---------------------------------------------------------------------------
// Angles between polylines
// ---------------------------------------------------------------------------

/**
 * Crossings of two polylines. Returns [{x, y, angle, i, j}] where angle (radians, 0..pi/2) is the
 * angle between the two (undirected) segments: pi/2 means perpendicular.
 */
export function polylineCrossings(a, b) {
    const out = [];
    for (let i = 0; i + 1 < a.length; i++) {
        for (let j = 0; j + 1 < b.length; j++) {
            const c = segmentCrossing(a[i], a[i + 1], b[j], b[j + 1]);
            if (c) out.push({ ...c, i, j });
        }
    }
    return out;
}

function segmentCrossing(p, q, a, b) {
    const rx = q[0] - p[0], ry = q[1] - p[1], sx = b[0] - a[0], sy = b[1] - a[1];
    const denom = rx * sy - ry * sx;
    if (denom === 0) return null;
    const apx = a[0] - p[0], apy = a[1] - p[1];
    const t = (apx * sy - apy * sx) / denom;
    const u = (apx * ry - apy * rx) / denom;
    if (t < 0 || t >= 1 || u < 0 || u >= 1) return null;
    // atan2 of |cross| and |dot| is accurate for nearly parallel as well as nearly perpendicular segments
    const angle = Math.atan2(Math.abs(denom), Math.abs(rx * sx + ry * sy));
    return { x: p[0] + t * rx, y: p[1] + t * ry, angle };
}

/**
 * The angle (radians, 0..pi/2, undirected) between two polylines at each of their crossings.
 * Convenience wrapper for the common case of one crossing: returns an array of angles.
 */
export function crossingAngles(a, b) {
    return polylineCrossings(a, b).map(c => c.angle);
}

/**
 * Crossings between two sets of polylines, restricted to segments touching `rect`
 * ({x0, y0, x1, y1}). Uses a uniform grid, so it is fast for hundreds of long lines.
 * Returns [{x, y, angle, a, b}] with a / b the indices into setA / setB.
 */
export function setCrossings(setA, setB, rect, cell = 20) {
    const grid = new Map();
    const key = (ix, iy) => ix * 100003 + iy;
    const touches = (p, q) =>
        Math.max(p[0], q[0]) >= rect.x0 && Math.min(p[0], q[0]) <= rect.x1 &&
        Math.max(p[1], q[1]) >= rect.y0 && Math.min(p[1], q[1]) <= rect.y1;
    const cellsOf = (p, q) => {
        const x0 = Math.max(Math.min(p[0], q[0]), rect.x0), x1 = Math.min(Math.max(p[0], q[0]), rect.x1);
        const y0 = Math.max(Math.min(p[1], q[1]), rect.y0), y1 = Math.min(Math.max(p[1], q[1]), rect.y1);
        const cells = [];
        for (let ix = Math.floor(x0 / cell); ix <= Math.floor(x1 / cell); ix++)
            for (let iy = Math.floor(y0 / cell); iy <= Math.floor(y1 / cell); iy++) cells.push(key(ix, iy));
        return cells;
    };
    setB.forEach((line, bi) => {
        const pts = line.points || line;
        for (let j = 0; j + 1 < pts.length; j++) {
            if (!touches(pts[j], pts[j + 1])) continue;
            for (const k of cellsOf(pts[j], pts[j + 1])) {
                if (!grid.has(k)) grid.set(k, []);
                grid.get(k).push([bi, j]);
            }
        }
    });
    const out = [];
    setA.forEach((line, ai) => {
        const pts = line.points || line;
        for (let i = 0; i + 1 < pts.length; i++) {
            if (!touches(pts[i], pts[i + 1])) continue;
            const seen = new Set();
            for (const k of cellsOf(pts[i], pts[i + 1])) {
                for (const [bi, j] of grid.get(k) || []) {
                    const id = bi * 1e7 + j;
                    if (seen.has(id)) continue;
                    seen.add(id);
                    const bpts = setB[bi].points || setB[bi];
                    const c = segmentCrossing(pts[i], pts[i + 1], bpts[j], bpts[j + 1]);
                    if (c) out.push({ ...c, a: ai, b: bi });
                }
            }
        }
    });
    return out;
}

// ---------------------------------------------------------------------------
// Tangency to a field
// ---------------------------------------------------------------------------

/**
 * Maximum misalignment (radians) between the direction of each segment of a polyline and the
 * field at the segment midpoint. `field(x, y)` returns {u, v}.
 *   oriented  - true: the polyline must run along the field (misalignment up to pi); false: a
 *               line that is only tangent, direction ignored (misalignment up to pi/2).
 *   sources   - segments whose midpoint is within `minDist` px of any source are skipped.
 *   minField  - segments where |field| < minField are skipped (near saddle points the field
 *               line direction changes too quickly to compare with a polyline chord).
 *   include   - optional predicate (x, y) => boolean selecting which midpoints are checked.
 * Returns {max, at: [x, y], n: number of segments checked}.
 */
export function maxTangentMisalignment(points, field, { sources = [], minDist = 5, minField = 0, oriented = true, include = null } = {}) {
    let max = 0, at = null, n = 0;
    for (let i = 0; i + 1 < points.length; i++) {
        const p = points[i], q = points[i + 1];
        const dx = q[0] - p[0], dy = q[1] - p[1];
        const len = Math.hypot(dx, dy);
        if (len === 0) continue;
        const mx = (p[0] + q[0]) / 2, my = (p[1] + q[1]) / 2;
        if (include && !include(mx, my)) continue;
        if (sources.some(s => Math.hypot(mx - s.x, my - s.y) < minDist)) continue;
        const f = field(mx, my);
        const fm = Math.hypot(f.u, f.v);
        if (fm === 0 || fm < minField) continue;
        const dot = dx * f.u + dy * f.v;
        const mis = Math.atan2(Math.abs(dx * f.v - dy * f.u), oriented ? dot : Math.abs(dot));
        n++;
        if (mis > max) { max = mis; at = [mx, my]; }
    }
    return { max, at, n };
}

// ---------------------------------------------------------------------------
// Comparing the same lines computed for different views
// ---------------------------------------------------------------------------

const inRect = (p, r) => p[0] >= r.x0 && p[0] <= r.x1 && p[1] >= r.y0 && p[1] <= r.y1;

/** Uniform-grid index of the segments of a polyline, for nearest-distance queries up to `cell` px. */
class SegmentIndex {
    constructor(segments, cell = 8) {
        this.cell = cell;
        this.grid = new Map();
        for (const [p, q] of segments) {
            const x0 = Math.floor((Math.min(p[0], q[0]) - cell) / cell), x1 = Math.floor((Math.max(p[0], q[0]) + cell) / cell);
            const y0 = Math.floor((Math.min(p[1], q[1]) - cell) / cell), y1 = Math.floor((Math.max(p[1], q[1]) + cell) / cell);
            // very long segments (rays to infinity) are kept in a separate short list
            if ((x1 - x0 + 1) * (y1 - y0 + 1) > 4000) { (this.long ||= []).push([p, q]); continue; }
            for (let ix = x0; ix <= x1; ix++) for (let iy = y0; iy <= y1; iy++) {
                const k = ix * 100003 + iy;
                if (!this.grid.has(k)) this.grid.set(k, []);
                this.grid.get(k).push([p, q]);
            }
        }
    }
    /** Distance from p to the polyline if below `cell`, else `cell` (a lower bound). */
    distance(p) {
        let best = this.cell;
        const k = Math.floor(p[0] / this.cell) * 100003 + Math.floor(p[1] / this.cell);
        for (const [a, b] of this.grid.get(k) || []) best = Math.min(best, pointSegmentDistance(p, a, b));
        for (const [a, b] of this.long || []) best = Math.min(best, pointSegmentDistance(p, a, b));
        return best;
    }
}

/**
 * Check that two computations of the same lines (e.g. with different `view`s) agree inside
 * `rect`. Lines are matched by (from, to, level); of the candidates the closest is taken. Each
 * line of one set with a vertex inside the rect must have a counterpart in the other set whose
 * in-rect part lies within `tol` px (vertices of each within `tol` of the other's polyline).
 * Distances above 8 px are reported as 8. Returns a list of problem strings.
 */
export function compareLinesInRect(linesA, linesB, rect, { tol = 0.5, levelTol = 1e-6 } = {}) {
    const problems = [];
    const sameLabel = (a, b) => a.from === b.from && a.to === b.to && Math.abs(a.level - b.level) <= levelTol;
    const describe = l => `[from ${l.from} to ${l.to}, level ${l.level}, ${l.points.length} pts, starts (${l.points[0].map(v => v.toFixed(1))})]`;
    const prep = lines => lines.map(l => {
        const inside = l.points.filter(p => inRect(p, rect));
        // only segments near the rect matter for distances from in-rect vertices
        const near = [];
        // a closed line is drawn with its closing segment (last point back to the first)
        const pts = l.closed && l.points.length > 2 ? [...l.points, l.points[0]] : l.points;
        for (let i = 0; i + 1 < pts.length; i++) {
            const p = pts[i], q = pts[i + 1];
            if (Math.max(p[0], q[0]) >= rect.x0 - 10 && Math.min(p[0], q[0]) <= rect.x1 + 10 &&
                Math.max(p[1], q[1]) >= rect.y0 - 10 && Math.min(p[1], q[1]) <= rect.y1 + 10) {
                near.push([p, q]);
            }
        }
        return { line: l, inside, index: new SegmentIndex(near) };
    });
    const A = prep(linesA), B = prep(linesB);
    const maxDist = (from, to) => from.inside.reduce((m, p) => Math.max(m, to.index.distance(p)), 0);
    const oneWay = (X, Y, nameX, nameY) => {
        for (const a of X) {
            if (a.inside.length === 0) continue;
            let best = Infinity, candidates = 0;
            for (const b of Y) {
                if (!sameLabel(a.line, b.line)) continue;
                candidates++;
                best = Math.min(best, Math.max(maxDist(a, b), maxDist(b, a)));
            }
            if (candidates === 0) problems.push(`${nameX} line ${describe(a.line)} has no ${nameY} line with the same from/to/level`);
            else if (best > tol) problems.push(`${nameX} line ${describe(a.line)} differs from its ${nameY} counterpart by ${best >= 8 ? '>= 8' : best.toFixed(2)} px inside the view (tolerance ${tol})`);
        }
    };
    oneWay(A, B, 'first', 'second');
    oneWay(B, A, 'second', 'first');
    return problems;
}

/**
 * Same idea for unlabelled contours (potential / magnetic): matched by level only and, among
 * equal levels, by geometry.
 */
export function compareContoursInRect(linesA, linesB, rect, { tol = 0.5, levelTol = 1e-6 } = {}) {
    const strip = l => ({ ...l, from: 0, to: 0 });
    return compareLinesInRect(linesA.map(strip), linesB.map(strip), rect, { tol, levelTol });
}
