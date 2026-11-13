// Field lines and contours for the analytic 2D field models in fieldmodels.js.
// See https://github.com/apontzen/maxwell/issues/3 (stage 2 and decisions D3-D5).
//
// These functions only compute geometry; draw.js does the drawing (clipping to the
// canvas, direction arrows, Gaussian-surface crossings, colours).
//
// CONTRACT
// --------
// All coordinates are canvas CSS pixels. `view` is {width, height} of the canvas.
// Returned polylines are NOT clipped to the view: lines are traced until they end on
// a charge or reach genuine infinity, independent of the view (rule: lines must not
// change configuration just because they leave the viewport). `view` may only be
// used to choose resolution (e.g. smaller steps where lines are visible).
//
// electricFieldLines(model, view) -> Array<{
//     points: Array<[x, y]>,   // ordered along the field direction (from + towards -)
//     from: number | 'infinity',  // index into model.sources, or 'infinity'
//     to:   number | 'infinity',
//     level: number,           // the value of the stream function psi on this line (mod q)
// }>
//   - model.kind === 'electrostatic_line_charge'.
//   - Lines are contours of model.streamFunctionAt at levels k / LINES_PER_UNIT_CHARGE (D4).
//     With the offset C = 0 (D5) the levels are exactly the multiples of 1/LINES_PER_UNIT_CHARGE.
//   - Every source of charge q has exactly LINES_PER_UNIT_CHARGE * |q| line ends attached, for q
//     a multiple of 1/LINES_PER_UNIT_CHARGE. (Lines from a charge back to the same charge are
//     impossible for a static field.)
//   - Includes lines from/to infinity, including those of the uniform field (model.uniform),
//     which may not touch any charge. Each physical line appears exactly once.
//   - A line ending on a charge ends exactly at the charge's position.
//
// magneticFieldLines(model, view) -> Array<{ points: Array<[x, y]>, closed: boolean, level: number }>
//   - model.kind === 'magnetostatic_wire'. Contours of model.potentialAt (= A_z) at levels
//     k * contourStep(model) (D4), covering at least the view; closed loops have closed = true
//     (first point repeated at the end is not required). Open contours (e.g. with a uniform field)
//     extend beyond the view.
//
// potentialContours(model, view) -> Array<{ points: Array<[x, y]>, closed: boolean, level: number }>
//   - model.kind === 'electrostatic_line_charge'. Contours of model.potentialAt at levels
//     k * contourStep(model), k != 0 plus the zero contour, covering at least the view.
//     draw.js colours level > 0 red, level < 0 blue, level == 0 grey (as the lecture text describes).
//
// contourStep(model) -> number
//   - The single level spacing used by both magneticFieldLines and potentialContours (D4), in
//     the units of model.potentialAt. Proportional to model.K so that the picture does not depend
//     on canvas size.
//
// IMPLEMENTATION (numerical core in levelsets2d.js)
// --------------
// All three are exact level sets, traced by predictor + Newton projection (no drift), with steps
// <= 2 px near the view and growing with the distance from it elsewhere; lines run out to
// "infinity" at 10 view diagonals. Stagnation points are found exactly (polynomial roots); a
// line/contour that runs exactly into one (by symmetry) is continued along the next branch, so
// each branch of the crossing is drawn once.
//   Electric: each source is seeded at the LPU*|q| exact crossings of psi = k/LPU on a small
//   circle around it. Lines from positive sources are traced along E and claim the slot of the
//   negative source they reach; the unclaimed slots of negative sources are traced backwards to
//   infinity. With a uniform field, lines touching no source are seeded where psi crosses a level
//   on the boundary of the view (plus SEED_MARGIN), skipping seeds on lines already traced.
//   Contours: seeded where the potential crosses a level on the boundary of the view (open
//   contours), and on the SOURCE_DISK_RADIUS circle and a ray from each source inside it (every
//   closed contour encloses a source); seeds on contours already traced are skipped. Contours are
//   traced through the disks around sources but the parts inside are dropped.
//   Coincident sources are traced as one source of their total charge (with zero-length lines
//   between the members' remaining line ends).


import { LINES_PER_UNIT_CHARGE } from './fieldmodels.js';
import {
    buildSystem, PsiEvaluator, PotentialEvaluator, stagnationPoints, softenStagnationPoints,
    prepareSnap, traceContour, sampleStep, seedBox, boxPoint, SegmentHash,
} from './levelsets2d.js';

export { LINES_PER_UNIT_CHARGE };

const LPU = LINES_PER_UNIT_CHARGE;
const TWO_PI = 2 * Math.PI;

// Contour spacing in units of model.K (see contourStep). A compromise: on a 668x335 canvas it
// gives 8 equipotentials around each charge of lecture 4 (+2,+2; 13 for +2,-2) and 3 closed
// field lines around each wire of the lecture 12 magnetic dipole (+-1, 80 px apart).
const CONTOUR_STEP_OVER_K = 0.5;
// Contours (magnetic field lines, equipotentials) stop at this distance from a source: the
// radius of the drawn charge / current circle (draw.js chargeSize)
export const SOURCE_DISK_RADIUS = 10;
// Seeds for lines/contours that touch no source are taken on the boundary of the view expanded by this
const SEED_MARGIN = 10;

// Diagnostics of the last call (for tests and debugging; not part of the contract)
export const diagnostics = {};

function wrap(a) {
    if (a > Math.PI) return a - TWO_PI * Math.ceil((a - Math.PI) / TWO_PI);
    if (a <= -Math.PI) return a + TWO_PI * Math.ceil((-Math.PI - a) / TWO_PI);
    return a;
}

export function contourStep(model) {
    return CONTOUR_STEP_OVER_K * model.K;
}

// ---------------------------------------------------------------------------------------------
// Electric field lines
// ---------------------------------------------------------------------------------------------

/** Radius around source i inside which psi is dominated by source i (monotonic on every circle). */
function dominanceRadius(T, i) {
    let s = Math.hypot(T.ux, T.uy) / T.K, dmin = Infinity;
    for (let j = 0; j < T.n; j++) {
        if (j === i) continue;
        const d = Math.hypot(T.xs[i] - T.xs[j], T.ys[i] - T.ys[j]);
        dmin = Math.min(dmin, d);
        s += Math.abs(T.qs[j]) / d;
    }
    let rho = 4.0;
    if (s > 0) rho = Math.min(rho, 0.35 * Math.abs(T.qs[i]) / s);
    rho = Math.min(rho, 0.3 * dmin);
    return Math.max(rho, 1e-9);
}

/**
 * The LPU*|q| points where the field lines of source i cross a small circle around it, with the
 * psi frame (unwrapped angles) at each.
 */
function chargeSlots(T, ev, i) {
    const q = T.qs[i];
    const nl = Math.round(LPU * Math.abs(q));
    if (nl === 0) return null;
    const rho = dominanceRadius(T, i);
    const xi = T.xs[i], yi = T.ys[i];
    const thc = new Float64Array(T.n);
    for (let j = 0; j < T.n; j++) if (j !== i) thc[j] = Math.atan2(yi - T.ys[j], xi - T.xs[j]);
    const c = 1 / (TWO_PI * T.K);
    const state = (alpha) => {
        const x = xi + rho * Math.cos(alpha), y = yi + rho * Math.sin(alpha);
        const th = new Float64Array(T.n);
        let s = 0;
        for (let j = 0; j < T.n; j++) {
            th[j] = j === i ? alpha : thc[j] + wrap(Math.atan2(y - T.ys[j], x - T.xs[j]) - thc[j]);
            s += T.qs[j] * th[j];
        }
        return { x, y, th, psi: s / TWO_PI + (T.ux * (y - T.myc) - T.uy * (x - T.mxc)) * c };
    };
    const a = state(-Math.PI).psi;
    const sgn = q > 0 ? 1 : -1;
    const k0 = q > 0 ? Math.ceil(LPU * a - 1e-9) : Math.floor(LPU * a + 1e-9);
    const list = [];
    for (let m = 0; m < nl; m++) {
        const level = (k0 + sgn * m) / LPU;
        // psi is monotonic in alpha (increasing for q > 0) on [-pi - 0.1, pi + 0.1]
        let lo = -Math.PI - 0.1, hi = Math.PI + 0.1;
        for (let it = 0; it < 60; it++) {
            const mid = 0.5 * (lo + hi);
            if ((state(mid).psi - level) * sgn < 0) lo = mid; else hi = mid;
        }
        const st = state(0.5 * (lo + hi));
        list.push({ angle: wrap(0.5 * (lo + hi)), x: st.x, y: st.y, th: st.th, level, line: -1 });
    }
    return { rho, list };
}

function circleCrossing(ax, ay, bx, by, cx, cy, r) {
    // parameter u in [0, 1] where segment a->b crosses the circle (a outside, b inside)
    const dx = bx - ax, dy = by - ay, fx = ax - cx, fy = ay - cy;
    const A = dx * dx + dy * dy, B = 2 * (fx * dx + fy * dy), C = fx * fx + fy * fy - r * r;
    const disc = Math.max(0, B * B - 4 * A * C);
    const u = A > 0 ? (-B - Math.sqrt(disc)) / (2 * A) : 0;
    return Math.min(1, Math.max(0, u));
}

export function electricFieldLines(model, view) {
    // Sources at the same point act as one source of the total charge. Trace that, then share
    // its line ends among the members and join the members' remaining ends with zero-length lines.
    const src = model.sources;
    const groups = [];
    const groupOf = new Int32Array(src.length).fill(-1);
    for (let i = 0; i < src.length; i++) {
        if (groupOf[i] >= 0) continue;
        groupOf[i] = groups.length;
        const g = [i];
        for (let j = i + 1; j < src.length; j++)
            if (groupOf[j] < 0 && Math.abs(src[j].x - src[i].x) < 1e-6 && Math.abs(src[j].y - src[i].y) < 1e-6) {
                groupOf[j] = groups.length; g.push(j);
            }
        groups.push(g);
    }
    if (groups.length === src.length) return electricLines(model, view);

    const merged = groups.map(g => ({ x: src[g[0]].x, y: src[g[0]].y, charge: g.reduce((a, i) => a + src[i].charge, 0) }));
    const keep = merged.map((c, k) => k).filter(k => merged[k].charge !== 0);
    const reduced = { ...model, sources: keep.map(k => merged[k]) };
    const lines = electricLines(reduced, view);
    // remaining capacity of each original source
    const cap = src.map(c => Math.round(LPU * Math.abs(c.charge)));
    const take = (g, sign) => {
        for (const i of groups[g]) if (cap[i] > 0 && Math.sign(src[i].charge) === sign) { cap[i]--; return i; }
        return -1;
    };
    for (const l of lines) {
        if (l.from !== 'infinity') l.from = take(keep[l.from], 1);
        if (l.to !== 'infinity') l.to = take(keep[l.to], -1);
    }
    for (let g = 0; g < groups.length; g++) {
        for (;;) {
            const a = groups[g].find(i => cap[i] > 0 && src[i].charge > 0);
            const b = groups[g].find(i => cap[i] > 0 && src[i].charge < 0);
            if (a === undefined || b === undefined) break;
            cap[a]--; cap[b]--;
            lines.push({ points: [[src[a].x, src[a].y], [src[b].x, src[b].y]], from: a, to: b, level: 0 });
        }
    }
    return lines;
}

function electricLines(model, view) {
    const T = buildSystem(model, view);
    const diag = { conflicts: 0, maxsteps: 0, freeHitCharge: 0, seeds: 0 };
    Object.keys(diagnostics).forEach(k => delete diagnostics[k]);
    Object.assign(diagnostics, diag);
    if (T.n === 0 && !T.hasUniform) return [];
    T.crit = stagnationPoints(T);
    const ev = new PsiEvaluator(T);
    const snaps = [];
    for (const s of T.crit) {
        const v = LPU * ev.principal(s.x, s.y);
        if (Math.abs(v - Math.round(v)) < 1e-8) {
            const p = prepareSnap(T, ev, s, 1 / LPU);
            if (p) snaps.push(p);
        }
    }
    T.snaps = snaps;

    const slots = [];
    for (let i = 0; i < T.n; i++) slots.push(chargeSlots(T, ev, i));
    const capturer = (sign) => (x, y) => {
        for (let j = 0; j < T.n; j++) {
            const sl = slots[j];
            if (!sl || T.qs[j] * sign <= 0) continue;
            const dx = x - T.xs[j], dy = y - T.ys[j];
            if (dx * dx + dy * dy < sl.rho * sl.rho) return j;
        }
        return -1;
    };
    const captureNeg = capturer(-1), capturePos = capturer(1);

    // Finish a line captured by source j: replace the final point by the slot's exact crossing of
    // the capture circle and the source position, and claim the slot.
    const arrive = (pts, j, lineIndex) => {
        const sl = slots[j];
        const [bx, by] = pts[pts.length - 1], [ax, ay] = pts[pts.length - 2];
        const u = circleCrossing(ax, ay, bx, by, T.xs[j], T.ys[j], sl.rho);
        const ang = Math.atan2(ay + u * (by - ay) - T.ys[j], ax + u * (bx - ax) - T.xs[j]);
        let best = 0, bestd = Infinity;
        sl.list.forEach((s, k) => {
            const d = Math.abs(wrap(ang - s.angle));
            if (d < bestd) { bestd = d; best = k; }
        });
        const slot = sl.list[best];
        if (slot.line >= 0) diag.conflicts++;
        slot.line = lineIndex;
        pts[pts.length - 1] = [slot.x, slot.y];
        pts.push([T.xs[j], T.ys[j]]);
    };

    const lines = [];
    // lines leaving positive charges, traced along E: to a negative charge or to infinity
    for (let i = 0; i < T.n; i++) {
        if (!slots[i] || T.qs[i] <= 0) continue;
        for (const slot of slots[i].list) {
            ev.copyStateFrom(slot.th);
            const res = traceContour(T, ev, slot.x, slot.y, slot.level, 1, { capture: captureNeg, snaps });
            const pts = [[T.xs[i], T.ys[i]], ...res.pts];
            const index = lines.length;
            slot.line = index;
            let to = 'infinity';
            if (res.end === 'capture') { arrive(pts, res.source, index); to = res.source; }
            else if (res.end !== 'infinity') diag.maxsteps++;
            lines.push({ points: pts, from: i, to, level: slot.level });
        }
    }
    // remaining lines of negative charges come from infinity: trace them backwards
    for (let i = 0; i < T.n; i++) {
        if (!slots[i] || T.qs[i] >= 0) continue;
        for (const slot of slots[i].list) {
            if (slot.line >= 0) continue;
            ev.copyStateFrom(slot.th);
            const res = traceContour(T, ev, slot.x, slot.y, slot.level, -1, { capture: capturePos, snaps });
            const pts = [[T.xs[i], T.ys[i]], ...res.pts];
            const index = lines.length;
            slot.line = index;
            let from = 'infinity';
            if (res.end === 'capture') { arrive(pts, res.source, index); from = res.source; diag.conflicts++; }
            else if (res.end !== 'infinity') diag.maxsteps++;
            pts.reverse();
            lines.push({ points: pts, from, to: i, level: slot.level });
        }
    }

    // lines of the uniform field that touch no charge: every one crossing the view crosses the
    // boundary of the (expanded) view, so seed on that boundary
    if (T.hasUniform) {
        const box = seedBox(T, SEED_MARGIN);
        const hash = new SegmentHash(box);
        for (const l of lines) hash.insertPolyline(l.points, 0);
        const seeds = psiBoxSeeds(T, ev, box);
        diag.seeds = seeds.length;
        for (const sd of seeds) {
            if (hash.contains(sd.x, sd.y)) continue;
            ev.copyStateFrom(sd.th);
            const f = traceContour(T, ev, sd.x, sd.y, sd.level, 1, { capture: captureNeg, snaps });
            if (f.end === 'capture') { diag.freeHitCharge++; continue; }
            ev.copyStateFrom(sd.th);
            const b = traceContour(T, ev, sd.x, sd.y, sd.level, -1, { capture: capturePos, snaps });
            if (b.end === 'capture') { diag.freeHitCharge++; continue; }
            if (f.end !== 'infinity' || b.end !== 'infinity') diag.maxsteps++;
            const pts = b.pts.reverse().concat(f.pts.slice(1));
            hash.insertPolyline(pts, 0);
            lines.push({ points: pts, from: 'infinity', to: 'infinity', level: sd.level });
        }
    }
    Object.assign(diagnostics, diag, { stagnation: T.crit.length, snaps: snaps.length });
    return lines;
}

/** Points on the box boundary where psi (unwrapped along the boundary) crosses a level k / LPU. */
function psiBoxSeeds(T, ev, box) {
    const P = 2 * ((box.x1 - box.x0) + (box.y1 - box.y0));
    let s = 0;
    let [x, y] = boxPoint(box, 0);
    ev.resetPrincipal(x, y);
    let v = ev.val(x, y);
    const seeds = [];
    while (s < P) {
        const thPrev = ev.saveState();
        const xp = x, yp = y, vp = v;
        const sp = s;
        s = Math.min(P, s + sampleStep(T, x, y));
        [x, y] = boxPoint(box, s);
        v = ev.val(x, y);
        let ka = Math.floor(LPU * vp), kb = Math.floor(LPU * v);
        if (ka === kb) continue;
        const lo = Math.min(ka, kb) + 1, hi = Math.max(ka, kb);
        for (let m = lo; m <= hi; m++) {
            const level = m / LPU;
            let t0 = sp, t1 = s;
            const f = (t) => { const [qx, qy] = boxPoint(box, t); return vp + ev.diff(qx, qy, xp, yp) - level; };
            const f0 = f(t0);
            for (let it = 0; it < 50; it++) {
                const tm = 0.5 * (t0 + t1);
                if ((f(tm) < 0) === (f0 < 0)) t0 = tm; else t1 = tm;
            }
            const [qx, qy] = boxPoint(box, 0.5 * (t0 + t1));
            const th = new Float64Array(T.n);
            for (let j = 0; j < T.n; j++)
                th[j] = thPrev[j] + wrap(Math.atan2(qy - T.ys[j], qx - T.xs[j]) - Math.atan2(yp - T.ys[j], xp - T.xs[j]));
            seeds.push({ x: qx, y: qy, th, level });
        }
    }
    return seeds;
}

// ---------------------------------------------------------------------------------------------
// Contours of the (single-valued) potential: magnetic field lines and equipotentials
// ---------------------------------------------------------------------------------------------

function seedsAlong(ev, step, pointAt, t0, t1, stepAt, out, extra) {
    // sample F along a curve, collecting crossings of the levels k * step
    let t = t0;
    let [x, y] = pointAt(t);
    let v = ev.val(x, y);
    while (t < t1) {
        const tp = t, vp = v;
        t = Math.min(t1, t + stepAt(t, x, y));
        [x, y] = pointAt(t);
        v = ev.val(x, y);
        const ka = Math.floor(vp / step), kb = Math.floor(v / step);
        if (ka === kb) continue;
        const lo = Math.min(ka, kb) + 1, hi = Math.max(ka, kb);
        for (let k = lo; k <= hi; k++) {
            const level = k * step;
            let a = tp, b = t, fa = vp - level;
            for (let it = 0; it < 50; it++) {
                const m = 0.5 * (a + b);
                const [mx, my] = pointAt(m);
                const fm = ev.val(mx, my) - level;
                if ((fm < 0) === (fa < 0)) { a = m; fa = fm; } else b = m;
            }
            const [sx, sy] = pointAt(0.5 * (a + b));
            out.push({ x: sx, y: sy, k, ...extra });
        }
    }
}

/** The point on the circle of radius r0 about source i where F = L, near (x, y) on that circle. */
function onCircle(T, ev, i, r0, x, y, L) {
    const cx = T.xs[i], cy = T.ys[i];
    let t = Math.atan2(y - cy, x - cx);
    const t0 = t;
    for (let it = 0; it < 8; it++) {
        const c = Math.cos(t), s = Math.sin(t);
        const r = ev.val(cx + r0 * c, cy + r0 * s) - L;
        const d = r0 * (c * ev.gy - s * ev.gx);
        if (!(Math.abs(d) > 0)) break;
        const dt = r / d;
        t -= Math.max(-0.05, Math.min(0.05, dt));
        if (Math.abs(dt) < 1e-13) break;
    }
    if (Math.abs(t - t0) > 0.1) t = t0;  // keep the chord estimate if Newton wandered
    return [cx + r0 * Math.cos(t), cy + r0 * Math.sin(t)];
}

function splitAtDisks(T, ev, L, path, closed, r0) {
    const n = path.length;
    const inside = new Int32Array(n).fill(-1);
    let any = false;
    for (let p = 0; p < n; p++) {
        const [x, y] = path[p];
        for (let i = 0; i < T.n; i++) {
            const dx = x - T.xs[i], dy = y - T.ys[i];
            if (dx * dx + dy * dy < r0 * r0) { inside[p] = i; any = true; break; }
        }
    }
    if (!any) return [{ points: path, closed }];
    let order = [];
    if (closed) {
        let start = 0;
        while (inside[start] < 0) start++;
        for (let k = 0; k <= n; k++) order.push((start + k) % n);
    } else {
        for (let k = 0; k < n; k++) order.push(k);
    }
    const pieces = [];
    let cur = null;
    for (let k = 0; k < order.length; k++) {
        const p = order[k];
        if (inside[p] < 0) {
            if (!cur) {
                cur = [];
                if (k > 0) {
                    const q = order[k - 1], i = inside[q];
                    const [bx, by] = path[q], [ax, ay] = path[p];   // from inside (b) to outside (a)
                    const u = circleCrossing(ax, ay, bx, by, T.xs[i], T.ys[i], r0);
                    cur.push(onCircle(T, ev, i, r0, ax + u * (bx - ax), ay + u * (by - ay), L));
                }
            }
            cur.push(path[p]);
        } else if (cur) {
            const q = order[k - 1], i = inside[p];
            const [ax, ay] = path[q], [bx, by] = path[p];
            const u = circleCrossing(ax, ay, bx, by, T.xs[i], T.ys[i], r0);
            cur.push(onCircle(T, ev, i, r0, ax + u * (bx - ax), ay + u * (by - ay), L));
            if (cur.length >= 2) pieces.push({ points: cur, closed: false });
            cur = null;
        }
    }
    if (cur && cur.length >= 2) pieces.push({ points: cur, closed: false });
    return pieces;
}

function potentialLevelSets(model, view, sigma) {
    const T = buildSystem(model, view);
    const step = contourStep(model);
    const r0 = SOURCE_DISK_RADIUS;
    const diag = { maxsteps: 0, seeds: 0, traced: 0 };
    T.crit = softenStagnationPoints(T, stagnationPoints(T));
    const ev = new PotentialEvaluator(T, sigma);
    const snapsByLevel = new Map();
    for (const s of T.crit) {
        const v = ev.val(s.x, s.y) / step;
        const k = Math.round(v);
        if (Math.abs(v - k) < 1e-9 * Math.max(1, Math.abs(v))) {
            const p = prepareSnap(T, ev, s, step);
            if (p) {
                if (!snapsByLevel.has(k)) snapsByLevel.set(k, []);
                snapsByLevel.get(k).push(p);
            }
        }
    }

    const box = seedBox(T, SEED_MARGIN);
    const seeds = [];
    // the boundary of the (expanded) view: every open contour crossing the view
    const P = 2 * ((box.x1 - box.x0) + (box.y1 - box.y0));
    seedsAlong(ev, step, (t) => boxPoint(box, t), 0, P, (t, x, y) => sampleStep(T, x, y), seeds);
    // every closed contour lying within the box encloses a source in the box, so crosses the
    // circle around it or a ray from it
    for (let i = 0; i < T.n; i++) {
        const xi = T.xs[i], yi = T.ys[i];
        if (xi < box.x0 || xi > box.x1 || yi < box.y0 || yi > box.y1) continue;
        seedsAlong(ev, step, (t) => [xi + r0 * Math.cos(t), yi + r0 * Math.sin(t)],
            -Math.PI, Math.PI, () => TWO_PI / 128, seeds);
        const dists = [xi - box.x0, box.x1 - xi, yi - box.y0, box.y1 - yi];
        const dirs = [[-1, 0], [1, 0], [0, -1], [0, 1]];
        let e = 0;
        for (let k = 1; k < 4; k++) if (dists[k] < dists[e]) e = k;
        if (dists[e] > r0) {
            const [dx, dy] = dirs[e];
            seedsAlong(ev, step, (t) => [xi + t * dx, yi + t * dy], r0, dists[e],
                (t, x, y) => sampleStep(T, x, y), seeds);
        }
    }
    diag.seeds = seeds.length;

    const insideDisk = (x, y) => {
        for (let i = 0; i < T.n; i++) {
            const dx = x - T.xs[i], dy = y - T.ys[i];
            if (dx * dx + dy * dy < r0 * r0 * (1 - 1e-9)) return true;
        }
        return false;
    };
    const hash = new SegmentHash(box);
    const out = [];
    for (const sd of seeds) {
        if (insideDisk(sd.x, sd.y)) continue;
        if (hash.contains(sd.x, sd.y, sd.k)) continue;
        diag.traced++;
        const L = sd.k * step;
        const snaps = snapsByLevel.get(sd.k) || [];
        const f = traceContour(T, ev, sd.x, sd.y, L, 1, { closeGate: true, snaps });
        let path, closed;
        if (f.end === 'closed') { path = f.pts; closed = true; }
        else {
            const b = traceContour(T, ev, sd.x, sd.y, L, -1, { snaps });
            if (f.end !== 'infinity' || b.end !== 'infinity') diag.maxsteps++;
            path = b.pts.reverse().concat(f.pts.slice(1));
            closed = false;
        }
        hash.insertPolyline(path, sd.k);
        if (closed) hash.insertPolyline([path[path.length - 1], path[0]], sd.k);
        // a level within rounding of 0 is reported as exactly 0
        const level = sd.k === 0 ? 0 : L;
        for (const piece of splitAtDisks(T, ev, L, path, closed, r0)) out.push({ ...piece, level });
    }
    Object.keys(diagnostics).forEach(k => delete diagnostics[k]);
    Object.assign(diagnostics, diag, { stagnation: T.crit.length });
    return out;
}

export function magneticFieldLines(model, view) {
    return potentialLevelSets(model, view, -1);
}

export function potentialContours(model, view) {
    return potentialLevelSets(model, view, 1);
}
