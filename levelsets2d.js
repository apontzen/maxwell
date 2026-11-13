// Numerical core for fieldlines2d.js: exact level-set tracing for the analytic 2D field models.
// See https://github.com/apontzen/maxwell/issues/3 (stage 2).
//
// Everything drawn in field-line mode is a level set of a harmonic function of the sources:
//   - electric field lines: contours of the stream function psi (multivalued: unwrapped here by
//     following each source's polar angle continuously along the path);
//   - magnetic field lines (A_z) and electric equipotentials (phi): contours of a single-valued
//     log potential.
// Contours are traced with a predictor (midpoint) step along the level set followed by a Newton
// projection back onto it, so points lie on the exact contour and there is no drift. Steps are
// small where the contour is visible and scale with the distance to the nearest source or
// stagnation point elsewhere, so tracing out to "infinity" is cheap.
//
// Stagnation points (zeros of the field) are found exactly as roots of a polynomial. A contour
// that runs exactly into a stagnation point (as happens by symmetry, e.g. on the axis between two
// equal charges) is snapped onto it and continued along the next outgoing branch, so that the
// branches of the crossing are each drawn once and the result does not depend on rounding.

import { SOFTEN } from './fieldmodels.js';

const TWO_PI = 2 * Math.PI;

export const VISIBLE_STEP = 2.0;      // px: maximum step where a contour is visible
const STEP_MARGIN = 24;               // px: region around the view in which VISIBLE_STEP applies
const FAR_STEP_FRACTION = 0.1;        // step grows as this fraction of the distance from that region
const NEAR_FRACTION = 0.1;            // step <= this fraction of the distance to a source / stagnation point
const MAX_TURN_COS = Math.cos(0.15);  // max turning angle of the tangent over a half step
const INFINITY_FACTOR = 10;           // "infinity" is this many view diagonals from the view centre
export const MAX_STEPS = 30000;

function wrap(a) {
    // wrap an angle difference into (-pi, pi]
    if (a > Math.PI) return a - TWO_PI * Math.ceil((a - Math.PI) / TWO_PI);
    if (a <= -Math.PI) return a + TWO_PI * Math.ceil((-Math.PI - a) / TWO_PI);
    return a;
}

/**
 * Collect everything the tracers need about a model and a view.
 */
export function buildSystem(model, view) {
    const src = model.sources;
    const n = src.length;
    const xs = new Float64Array(n), ys = new Float64Array(n), qs = new Float64Array(n);
    for (let i = 0; i < n; i++) { xs[i] = src[i].x; ys[i] = src[i].y; qs[i] = src[i].charge; }
    const vw = (view && view.width) || model.width, vh = (view && view.height) || model.height;
    const diag = Math.hypot(vw, vh);
    const T = {
        model, n, xs, ys, qs,
        K: model.K, ux: model.uniform.x, uy: model.uniform.y,
        // the centre and length used by the model's normalisation of psi and phi
        mxc: model.width / 2, myc: model.height / 2, lnW2: Math.log(model.width * model.width),
        vw, vh, diag,
        cx: vw / 2, cy: vh / 2,
        Rinf: INFINITY_FACTOR * diag,
        sx0: -STEP_MARGIN, sy0: -STEP_MARGIN, sx1: vw + STEP_MARGIN, sy1: vh + STEP_MARGIN,
        crit: [],          // stagnation points used to limit the step
        snaps: [],         // stagnation points lying exactly on a level traced in this run
        hasUniform: model.uniform.x !== 0 || model.uniform.y !== 0,
    };
    return T;
}

// ---------------------------------------------------------------------------------------------
// Evaluators. Each provides val(x, y) -> F (also setting gx, gy = grad F) and grad(x, y) (sets
// gx, gy only). The tangent followed in the "forward" direction is rot_-90(grad F) = (gy, -gx).
// ---------------------------------------------------------------------------------------------

/** Stream function psi (electric field lines); forward tangent = E (unsoftened). */
export class PsiEvaluator {
    constructor(T) {
        this.T = T;
        this.th = new Float64Array(T.n);  // unwrapped polar angle of the point about each source
        this.pa = new Float64Array(T.n);  // the principal value at the last evaluation
        this.gx = 0; this.gy = 0;
        this.c = 1 / (TWO_PI * T.K);
    }
    /** Set the unwrapped angles to their principal values at (x, y) (the frame of model.streamFunctionAt). */
    resetPrincipal(x, y) {
        const T = this.T;
        for (let i = 0; i < T.n; i++) {
            const a = Math.atan2(y - T.ys[i], x - T.xs[i]);
            this.th[i] = a; this.pa[i] = a;
        }
    }
    copyStateFrom(th) {
        const T = this.T;
        for (let i = 0; i < T.n; i++) { this.th[i] = th[i]; this.pa[i] = wrap(th[i]); }
    }
    saveState() { return Float64Array.from(this.th); }
    grad(x, y) {
        const T = this.T, xs = T.xs, ys = T.ys, qs = T.qs;
        let ex = 0, ey = 0;
        for (let i = 0; i < T.n; i++) {
            const dx = x - xs[i], dy = y - ys[i];
            const f = qs[i] / (dx * dx + dy * dy);
            ex += f * dx; ey += f * dy;
        }
        // E / K = sum q r/r^2 + u/K ; grad psi = rot90(E) / (2 pi K)
        ex = ex / TWO_PI + T.ux * this.c;
        ey = ey / TWO_PI + T.uy * this.c;
        this.gx = -ey; this.gy = ex;
    }
    val(x, y) {
        const T = this.T, xs = T.xs, ys = T.ys, qs = T.qs, th = this.th, pa = this.pa;
        let ex = 0, ey = 0, s = 0;
        for (let i = 0; i < T.n; i++) {
            const dx = x - xs[i], dy = y - ys[i];
            const q = qs[i];
            const f = q / (dx * dx + dy * dy);
            ex += f * dx; ey += f * dy;
            const a = Math.atan2(dy, dx);
            let d = a - pa[i];
            if (d > Math.PI) d -= TWO_PI; else if (d < -Math.PI) d += TWO_PI;
            th[i] += d; pa[i] = a;
            s += q * th[i];
        }
        ex = ex / TWO_PI + T.ux * this.c;
        ey = ey / TWO_PI + T.uy * this.c;
        this.gx = -ey; this.gy = ex;
        return s / TWO_PI + (T.ux * (y - T.myc) - T.uy * (x - T.mxc)) * this.c;
    }
    /** psi(x, y) - psi(x0, y0) along a short path (no branch cut crossed). */
    diff(x, y, x0, y0) {
        const T = this.T;
        let s = 0;
        for (let i = 0; i < T.n; i++)
            s += T.qs[i] * wrap(Math.atan2(y - T.ys[i], x - T.xs[i]) - Math.atan2(y0 - T.ys[i], x0 - T.xs[i]));
        return s / TWO_PI + (T.ux * (y - y0) - T.uy * (x - x0)) * this.c;
    }
    principal(x, y) {
        const T = this.T;
        let s = 0;
        for (let i = 0; i < T.n; i++) s += T.qs[i] * Math.atan2(y - T.ys[i], x - T.xs[i]);
        return s / TWO_PI + (T.ux * (y - T.myc) - T.uy * (x - T.mxc)) * this.c;
    }
}

/**
 * The model's (softened) potential: sigma = +1 gives phi (electric model.potentialAt),
 * sigma = -1 gives A_z (magnetic model.potentialAt). grad F = -sigma * E_equivalent (softened),
 * so for A_z the forward tangent (gy, -gx) is B = model.fieldAt.
 */
export class PotentialEvaluator {
    constructor(T, sigma) {
        this.T = T; this.sigma = sigma;
        this.gx = 0; this.gy = 0;
        this.halfK = 0.5 * T.K;
    }
    grad(x, y) {
        const T = this.T, xs = T.xs, ys = T.ys, qs = T.qs;
        let ex = 0, ey = 0;
        for (let i = 0; i < T.n; i++) {
            const dx = x - xs[i], dy = y - ys[i];
            const f = qs[i] / (dx * dx + dy * dy + SOFTEN);
            ex += f * dx; ey += f * dy;
        }
        ex = T.K * ex + T.ux; ey = T.K * ey + T.uy;
        this.gx = -this.sigma * ex; this.gy = -this.sigma * ey;
    }
    val(x, y) {
        const T = this.T, xs = T.xs, ys = T.ys, qs = T.qs;
        let ex = 0, ey = 0, l = 0;
        for (let i = 0; i < T.n; i++) {
            const dx = x - xs[i], dy = y - ys[i];
            const r2 = dx * dx + dy * dy + SOFTEN;
            const f = qs[i] / r2;
            ex += f * dx; ey += f * dy;
            l += qs[i] * (Math.log(r2) - T.lnW2);
        }
        ex = T.K * ex + T.ux; ey = T.K * ey + T.uy;
        this.gx = -this.sigma * ex; this.gy = -this.sigma * ey;
        const phi = -(T.ux * (x - T.mxc) + T.uy * (y - T.myc)) - this.halfK * l;
        return this.sigma * phi;
    }
    diff(x, y, x0, y0) { return this.val(x, y) - this.val(x0, y0); }
}

// ---------------------------------------------------------------------------------------------
// Stagnation points: zeros of E = K sum q_i / conj(z - z_i) + U, i.e. roots of
// P(z) = conj(U)/K prod (z - z_j) + sum_i q_i prod_{j != i} (z - z_j).
// ---------------------------------------------------------------------------------------------

function polyMulLinear(re, im, ar, ai) {
    // (sum c_k z^k) * (z - a)
    const n = re.length;
    const r = new Float64Array(n + 1), i = new Float64Array(n + 1);
    for (let k = 0; k < n; k++) {
        r[k + 1] += re[k]; i[k + 1] += im[k];
        r[k] -= re[k] * ar - im[k] * ai;
        i[k] -= re[k] * ai + im[k] * ar;
    }
    return [r, i];
}

function polyRoots(re, im) {
    // Durand-Kerner on a polynomial with coefficients a_0..a_d (a_d != 0)
    const d = re.length - 1;
    if (d < 1) return [];
    const lr = re[d], li = im[d], l2 = lr * lr + li * li;
    const ar = new Float64Array(d + 1), ai = new Float64Array(d + 1);
    let bound = 0;
    for (let k = 0; k <= d; k++) {
        ar[k] = (re[k] * lr + im[k] * li) / l2;
        ai[k] = (im[k] * lr - re[k] * li) / l2;
        if (k < d) bound = Math.max(bound, Math.hypot(ar[k], ai[k]));
    }
    const R = 1 + bound;
    const zr = new Float64Array(d), zi = new Float64Array(d);
    for (let k = 0; k < d; k++) {
        const t = TWO_PI * k / d + 0.4;
        zr[k] = R * 0.9 * Math.cos(t); zi[k] = R * 0.9 * Math.sin(t);
    }
    for (let it = 0; it < 1000; it++) {
        let change = 0;
        for (let k = 0; k < d; k++) {
            // p(z_k) by Horner
            let pr = 1, pi = 0;
            for (let j = d - 1; j >= 0; j--) {
                const nr = pr * zr[k] - pi * zi[k] + ar[j];
                const ni = pr * zi[k] + pi * zr[k] + ai[j];
                pr = nr; pi = ni;
            }
            let qr = 1, qi = 0;
            for (let j = 0; j < d; j++) {
                if (j === k) continue;
                const dr = zr[k] - zr[j], di = zi[k] - zi[j];
                const nr = qr * dr - qi * di, ni = qr * di + qi * dr;
                qr = nr; qi = ni;
            }
            const q2 = qr * qr + qi * qi;
            if (q2 === 0) { zr[k] += 1e-8; continue; }
            const wr = (pr * qr + pi * qi) / q2, wi = (pi * qr - pr * qi) / q2;
            zr[k] -= wr; zi[k] -= wi;
            change = Math.max(change, Math.hypot(wr, wi) / (1 + Math.hypot(zr[k], zi[k])));
        }
        if (change < 1e-15) break;
    }
    const out = [];
    for (let k = 0; k < d; k++) out.push([zr[k], zi[k]]);
    return out;
}

/**
 * Stagnation points of the unsoftened field, as [{x, y, m}] with m the multiplicity.
 */
export function stagnationPoints(T) {
    // sources with zero charge don't contribute (and would add a spurious root at their position)
    const idx = [];
    for (let j = 0; j < T.n; j++) if (T.qs[j] !== 0) idx.push(j);
    const n = idx.length;
    if (n === 0) return [];
    const L = T.diag, cx = T.cx, cy = T.cy;
    const zr = new Float64Array(n), zi = new Float64Array(n), qs = new Float64Array(n);
    for (let j = 0; j < n; j++) { zr[j] = (T.xs[idx[j]] - cx) / L; zi[j] = (T.ys[idx[j]] - cy) / L; qs[j] = T.qs[idx[j]]; }
    const cr = T.ux * L / T.K, ci = -T.uy * L / T.K;   // conj(U) L / K
    let all = [new Float64Array([1]), new Float64Array([0])];
    for (let j = 0; j < n; j++) all = polyMulLinear(all[0], all[1], zr[j], zi[j]);
    const Pr = new Float64Array(n + 1), Pi = new Float64Array(n + 1);
    for (let k = 0; k <= n; k++) {
        Pr[k] = all[0][k] * cr - all[1][k] * ci;
        Pi[k] = all[0][k] * ci + all[1][k] * cr;
    }
    for (let i = 0; i < n; i++) {
        let p = [new Float64Array([1]), new Float64Array([0])];
        for (let j = 0; j < n; j++) if (j !== i) p = polyMulLinear(p[0], p[1], zr[j], zi[j]);
        for (let k = 0; k < n; k++) { Pr[k] += qs[i] * p[0][k]; Pi[k] += qs[i] * p[1][k]; }
    }
    let maxc = 0;
    for (let k = 0; k <= n; k++) maxc = Math.max(maxc, Math.hypot(Pr[k], Pi[k]));
    if (maxc === 0) return [];
    let d = n;
    while (d > 0 && Math.hypot(Pr[d], Pi[d]) < 1e-13 * maxc) d--;
    if (d < 1) return [];
    const roots = polyRoots(Pr.slice(0, d + 1), Pi.slice(0, d + 1));

    // polish on g(z) = sum q / (z - z_j) + c, then cluster (multiple roots)
    const g = (xr, xi) => {
        let gr = cr, gi = ci, dr = 0, di = 0;
        for (let j = 0; j < n; j++) {
            const ar = xr - zr[j], ai = xi - zi[j];
            const a2 = ar * ar + ai * ai;
            const ir = ar / a2, ii = -ai / a2;          // 1 / (z - z_j)
            gr += qs[j] * ir; gi += qs[j] * ii;
            dr -= qs[j] * (ir * ir - ii * ii); di -= qs[j] * 2 * ir * ii;   // -q/(z-z_j)^2
        }
        return [gr, gi, dr, di];
    };
    const found = [];
    for (const [r0, i0] of roots) {
        if (!isFinite(r0) || !isFinite(i0) || Math.hypot(r0, i0) > 4 * INFINITY_FACTOR) continue;
        found.push({ zr: r0, zi: i0 });
    }
    // cluster
    const used = new Array(found.length).fill(false);
    const out = [];
    for (let a = 0; a < found.length; a++) {
        if (used[a]) continue;
        const members = [found[a]]; used[a] = true;
        for (let b = a + 1; b < found.length; b++) {
            if (used[b]) continue;
            if (Math.hypot(found[b].zr - found[a].zr, found[b].zi - found[a].zi) < 3e-6) {
                members.push(found[b]); used[b] = true;
            }
        }
        let xr = 0, xi = 0;
        for (const m of members) { xr += m.zr; xi += m.zi; }
        xr /= members.length; xi /= members.length;
        if (members.length === 1) {
            for (let it = 0; it < 30; it++) {
                const [gr, gi, dr, di] = g(xr, xi);
                const d2 = dr * dr + di * di;
                if (d2 === 0) break;
                const sr = (gr * dr + gi * di) / d2, si = (gi * dr - gr * di) / d2;
                xr -= sr; xi -= si;
                if (Math.hypot(sr, si) < 1e-16) break;
            }
        }
        out.push({ x: cx + L * xr, y: cy + L * xi, m: members.length });
    }
    // drop anything numerically on top of a source (cannot happen for distinct sources)
    return out.filter(s => {
        if (!isFinite(s.x) || !isFinite(s.y)) return false;
        for (let j = 0; j < T.n; j++) if (Math.hypot(s.x - T.xs[j], s.y - T.ys[j]) < 1e-9 * L) return false;
        return true;
    });
}

/** Move unsoftened stagnation points to the zeros of the softened field (2D Newton). */
export function softenStagnationPoints(T, pts) {
    return pts.map(p => {
        if (p.m !== 1) return { ...p };
        let x = p.x, y = p.y;
        for (let it = 0; it < 30; it++) {
            let ex = T.ux, ey = T.uy, jxx = 0, jxy = 0, jyy = 0;
            for (let i = 0; i < T.n; i++) {
                const dx = x - T.xs[i], dy = y - T.ys[i];
                const r2 = dx * dx + dy * dy + SOFTEN;
                const kq = T.K * T.qs[i];
                ex += kq * dx / r2; ey += kq * dy / r2;
                const r4 = r2 * r2;
                jxx += kq * (r2 - 2 * dx * dx) / r4;
                jxy += kq * (-2 * dx * dy) / r4;
                jyy += kq * (r2 - 2 * dy * dy) / r4;
            }
            const det = jxx * jyy - jxy * jxy;
            if (det === 0) break;
            const sx = (jyy * ex - jxy * ey) / det, sy = (jxx * ey - jxy * ex) / det;
            x -= sx; y -= sy;
            if (Math.hypot(sx, sy) < 1e-12 * T.diag) break;
        }
        if (!isFinite(x) || !isFinite(y) || Math.hypot(x - p.x, y - p.y) > 2) return { ...p };
        return { x, y, m: 1 };
    });
}

/**
 * Prepare a stagnation point that lies exactly on a traced level for snapping: its snapping
 * radius and the angles of its branches, sorted, each flagged incoming/outgoing for the forward
 * direction. Returns null if the branch structure can't be resolved (then it is not snapped).
 */
export function prepareSnap(T, ev, s, levelSpacing) {
    // distance to the nearest source or other stagnation point
    let ds = Infinity;
    for (let i = 0; i < T.n; i++) ds = Math.min(ds, Math.hypot(s.x - T.xs[i], s.y - T.ys[i]));
    for (const o of T.crit) if (o !== s) ds = Math.min(ds, Math.hypot(s.x - o.x, s.y - o.y) || Infinity);
    // curvature of F at the saddle: keep the variation of F over the snap radius << level spacing
    const h = 1e-3 * Math.max(1e-6, Math.min(ds, 10));
    ev.grad(s.x + h, s.y); const ax = ev.gx, ay = ev.gy;
    ev.grad(s.x - h, s.y); const bx = ev.gx, by = ev.gy;
    const hess = Math.hypot(ax - bx, ay - by) / (2 * h);
    let rs = Math.min(1.0, 0.05 * ds);
    if (hess > 0) rs = Math.min(rs, Math.sqrt(0.01 * levelSpacing / hess));
    const nb = 2 * (s.m + 1);
    const ns = 24 * nb;
    const vals = new Float64Array(ns + 1);
    const pt = (t) => [s.x + rs * Math.cos(t), s.y + rs * Math.sin(t)];
    for (let k = 0; k <= ns; k++) {
        const [x, y] = pt(TWO_PI * k / ns - Math.PI);
        vals[k] = ev.diff(x, y, s.x, s.y);
    }
    const branches = [];
    for (let k = 0; k < ns; k++) {
        const a = vals[k], b = vals[k + 1];
        if ((a < 0) === (b < 0) && a !== 0) continue;
        if (b === 0) continue; // counted in the next interval
        let t0 = TWO_PI * k / ns - Math.PI, t1 = TWO_PI * (k + 1) / ns - Math.PI, f0 = a;
        for (let it = 0; it < 50; it++) {
            const tm = 0.5 * (t0 + t1);
            const [x, y] = pt(tm);
            const fm = ev.diff(x, y, s.x, s.y);
            if ((fm < 0) === (f0 < 0)) { t0 = tm; f0 = fm; } else t1 = tm;
        }
        const t = 0.5 * (t0 + t1);
        const [x, y] = pt(t);
        ev.grad(x, y);
        const out = (ev.gy * Math.cos(t) - ev.gx * Math.sin(t)) > 0;  // forward tangent . radial
        branches.push({ angle: t, out });
    }
    if (branches.length !== nb) return null;
    for (let k = 0; k < nb; k++) if (branches[k].out === branches[(k + 1) % nb].out) return null;
    return { x: s.x, y: s.y, rs, branches };
}

// ---------------------------------------------------------------------------------------------
// The tracer
// ---------------------------------------------------------------------------------------------

function stepLimit(T, x, y) {
    const bx = Math.max(T.sx0 - x, 0, x - T.sx1), by = Math.max(T.sy0 - y, 0, y - T.sy1);
    let h = VISIBLE_STEP + FAR_STEP_FRACTION * Math.sqrt(bx * bx + by * by);
    let d2 = Infinity;
    for (let i = 0; i < T.n; i++) {
        const dx = x - T.xs[i], dy = y - T.ys[i];
        const r2 = dx * dx + dy * dy;
        if (r2 < d2) d2 = r2;
    }
    const crit = T.crit;
    for (let i = 0; i < crit.length; i++) {
        const dx = x - crit[i].x, dy = y - crit[i].y;
        const r2 = dx * dx + dy * dy;
        if (r2 < d2) d2 = r2;
    }
    const hn = NEAR_FRACTION * Math.sqrt(d2);
    if (hn < h) h = hn;
    return h > 1e-9 ? h : 1e-9;
}

/**
 * Trace the contour F = L from (x, y) (on the contour; the evaluator state must be set for it) in
 * direction dir (+1 forward along (gy, -gx), -1 backward).
 *
 * opts:
 *   capture(x, y) -> index of a source whose capture disk contains (x, y), or -1  (electric)
 *   closeGate: true to stop when the contour returns to its start (closed loops)
 *   snaps: prepared stagnation points to snap onto
 *
 * Returns {pts: [[x,y],...] (including the start), end: 'capture'|'infinity'|'closed'|'maxsteps',
 *          source (for capture)}
 */
export function traceContour(T, ev, x, y, L, dir, opts) {
    const pts = [[x, y]];
    const snaps = opts.snaps || [];
    const capture = opts.capture;
    const closeGate = opts.closeGate;
    const R2 = T.Rinf * T.Rinf;
    ev.val(x, y);
    let gx = ev.gx, gy = ev.gy;
    const x0 = x, y0 = y;
    let t0x = 0, t0y = 0;
    {
        const gn = Math.hypot(gx, gy);
        if (gn > 0) { t0x = dir * gy / gn; t0y = -dir * gx / gn; }
    }
    let arc = 0;
    let leaving = null;
    for (let step = 0; step < MAX_STEPS; step++) {
        let gn = Math.hypot(gx, gy);
        if (!(gn > 0)) return { pts, end: 'maxsteps' };
        const tx = dir * gy / gn, ty = -dir * gx / gn;
        let h = stepLimit(T, x, y);
        let mx = 0, my = 0;
        for (let tries = 0; tries < 30; tries++) {
            ev.grad(x + 0.5 * h * tx, y + 0.5 * h * ty);
            const mn = Math.hypot(ev.gx, ev.gy);
            mx = dir * ev.gy / mn; my = -dir * ev.gx / mn;
            if (tx * mx + ty * my >= MAX_TURN_COS || h < 1e-9) break;
            h *= 0.5;
        }
        let xn = x + h * mx, yn = y + h * my;
        for (let it = 0; it < 6; it++) {
            const r = ev.val(xn, yn) - L;
            const g2 = ev.gx * ev.gx + ev.gy * ev.gy;
            if (!(g2 > 0)) break;
            const c = r / g2;
            xn -= c * ev.gx; yn -= c * ev.gy;
            if (Math.abs(r) < 1e-7 * h * Math.sqrt(g2)) break;
        }
        gx = ev.gx; gy = ev.gy;
        const px = x, py = y;
        x = xn; y = yn;
        arc += h;
        pts.push([x, y]);

        if (capture) {
            const j = capture(x, y);
            if (j >= 0) return { pts, end: 'capture', source: j };
        }

        if (leaving) {
            if ((x - leaving.x) ** 2 + (y - leaving.y) ** 2 > 2.25 * leaving.rs * leaving.rs) leaving = null;
        }
        for (let k = 0; k < snaps.length; k++) {
            const s = snaps[k];
            if (s === leaving) continue;
            const dx = x - s.x, dy = y - s.y;
            if (dx * dx + dy * dy >= s.rs * s.rs) continue;
            // arrived at a stagnation point on this contour: continue along the paired branch
            const beta = Math.atan2(dy, dx);
            const want = dir > 0 ? false : true; // arriving forward along an incoming branch
            const nb = s.branches.length;
            let best = -1, bestd = Infinity;
            for (let b = 0; b < nb; b++) {
                if (s.branches[b].out !== want) continue;
                const d = Math.abs(wrap(beta - s.branches[b].angle));
                if (d < bestd) { bestd = d; best = b; }
            }
            if (best < 0) continue;
            const exit = s.branches[(best + (dir > 0 ? 1 : nb - 1)) % nb];
            pts[pts.length - 1] = [s.x, s.y];
            ev.val(s.x, s.y);
            x = s.x + s.rs * Math.cos(exit.angle);
            y = s.y + s.rs * Math.sin(exit.angle);
            for (let it = 0; it < 6; it++) {
                const r = ev.val(x, y) - L;
                const g2 = ev.gx * ev.gx + ev.gy * ev.gy;
                if (!(g2 > 0)) break;
                x -= r * ev.gx / g2; y -= r * ev.gy / g2;
                if (Math.abs(r) < 1e-9 * s.rs * Math.sqrt(g2)) break;
            }
            ev.val(x, y);
            gx = ev.gx; gy = ev.gy;
            pts.push([x, y]);
            leaving = s;
            break;
        }

        if (closeGate && arc > 4 * h) {
            const sp = (px - x0) * t0x + (py - y0) * t0y, sn = (x - x0) * t0x + (y - y0) * t0y;
            if (sp < 0 && sn >= 0) {
                const d2 = (x - x0) ** 2 + (y - y0) ** 2;
                if (d2 < 9 * h * h) {
                    pts.pop();
                    return { pts, end: 'closed' };
                }
            }
        }

        const rx = x - T.cx, ry = y - T.cy;
        if (rx * rx + ry * ry > R2 && rx * gy * dir - ry * gx * dir > 0) return { pts, end: 'infinity' };
    }
    return { pts, end: 'maxsteps' };
}

// ---------------------------------------------------------------------------------------------
// Seeding helpers
// ---------------------------------------------------------------------------------------------

/** Step for sampling a seed curve at (x, y). */
export function sampleStep(T, x, y) {
    let h = VISIBLE_STEP;
    let d2 = Infinity;
    for (let i = 0; i < T.n; i++) {
        const dx = x - T.xs[i], dy = y - T.ys[i];
        d2 = Math.min(d2, dx * dx + dy * dy);
    }
    for (const c of T.crit) d2 = Math.min(d2, (x - c.x) ** 2 + (y - c.y) ** 2);
    h = Math.min(h, NEAR_FRACTION * Math.sqrt(d2));
    return Math.max(h, 1e-4);
}

/**
 * Choose a seed box: the view expanded by about `margin`, nudged so that no source lies within
 * 0.5 px of its edge.
 */
export function seedBox(T, margin) {
    for (let k = 0; k < 20; k++) {
        const m = margin + 1.37 * k;
        const b = { x0: -m, y0: -m, x1: T.vw + m, y1: T.vh + m };
        let ok = true;
        for (let i = 0; i < T.n && ok; i++) {
            const x = T.xs[i], y = T.ys[i];
            const nearX = (Math.abs(x - b.x0) < 0.5 || Math.abs(x - b.x1) < 0.5) && y > b.y0 - 0.5 && y < b.y1 + 0.5;
            const nearY = (Math.abs(y - b.y0) < 0.5 || Math.abs(y - b.y1) < 0.5) && x > b.x0 - 0.5 && x < b.x1 + 0.5;
            if (nearX || nearY) ok = false;
        }
        if (ok) return b;
    }
    return { x0: -margin, y0: -margin, x1: T.vw + margin, y1: T.vh + margin };
}

/** Point at perimeter parameter s along the box (clockwise on screen from the top-left corner). */
export function boxPoint(b, s) {
    const w = b.x1 - b.x0, h = b.y1 - b.y0;
    if (s < w) return [b.x0 + s, b.y0];
    s -= w;
    if (s < h) return [b.x1, b.y0 + s];
    s -= h;
    if (s < w) return [b.x1 - s, b.y1];
    s -= w;
    return [b.x0, b.y1 - Math.min(s, h)];
}

// ---------------------------------------------------------------------------------------------
// Spatial hash of traced segments, used to recognise seeds lying on contours already traced.
// ---------------------------------------------------------------------------------------------

export class SegmentHash {
    constructor(b, cell = 8) {
        this.x0 = b.x0 - 2 * cell; this.y0 = b.y0 - 2 * cell;
        this.cell = cell;
        this.nx = Math.ceil((b.x1 - b.x0) / cell) + 4;
        this.ny = Math.ceil((b.y1 - b.y0) / cell) + 4;
        this.cells = new Map();
    }
    insertPolyline(pts, tag) {
        const c = this.cell;
        for (let k = 0; k + 1 < pts.length; k++) {
            const [ax, ay] = pts[k], [bx, by] = pts[k + 1];
            let i0 = Math.floor((Math.min(ax, bx) - this.x0) / c), i1 = Math.floor((Math.max(ax, bx) - this.x0) / c);
            let j0 = Math.floor((Math.min(ay, by) - this.y0) / c), j1 = Math.floor((Math.max(ay, by) - this.y0) / c);
            if (i1 < 0 || j1 < 0 || i0 >= this.nx || j0 >= this.ny) continue;
            i0 = Math.max(i0, 0); j0 = Math.max(j0, 0);
            i1 = Math.min(i1, this.nx - 1); j1 = Math.min(j1, this.ny - 1);
            if ((i1 - i0 + 1) * (j1 - j0 + 1) > 64) continue; // long segment far from any seed: skip
            const seg = [ax, ay, bx, by, tag];
            for (let i = i0; i <= i1; i++) for (let j = j0; j <= j1; j++) {
                const key = i * this.ny + j;
                let list = this.cells.get(key);
                if (!list) { list = []; this.cells.set(key, list); }
                list.push(seg);
            }
        }
    }
    /** Is (x, y) on a stored segment (with the given tag, if tag !== undefined)? */
    contains(x, y, tag) {
        const i = Math.floor((x - this.x0) / this.cell), j = Math.floor((y - this.y0) / this.cell);
        if (i < 0 || j < 0 || i >= this.nx || j >= this.ny) return false;
        const list = this.cells.get(i * this.ny + j);
        if (!list) return false;
        for (const [ax, ay, bx, by, t] of list) {
            if (tag !== undefined && t !== tag) continue;
            const dx = bx - ax, dy = by - ay;
            const l2 = dx * dx + dy * dy;
            let u = l2 > 0 ? ((x - ax) * dx + (y - ay) * dy) / l2 : 0;
            if (u < -1e-9 || u > 1 + 1e-9) continue;
            u = Math.min(1, Math.max(0, u));
            const ex = ax + u * dx - x, ey = ay + u * dy - y;
            // tolerance: the chord error of a step that turns by < ~0.3 rad, but much less than the
            // distance between distinct contours of the same level
            const tol = 0.05 * Math.sqrt(l2) + 1e-7;
            if (ex * ex + ey * ey <= tol * tol) return true;
        }
        return false;
    }
}
