// The field-line rules, as functions returning a list of human-readable problem strings (an empty
// list means the rule holds). Environment-independent: the caller passes in the modules under
// test, so the same code runs in the browser (fieldline-rules.spec.js imports /fieldmodels.js and
// /fieldlines2d.js in the page and calls these) and in Node (fieldline-checks.spec.js runs them
// against an exact reference implementation to prove the rules are satisfiable).
//
//   env = { models: <fieldmodels.js>, lines2d: <fieldlines2d.js>, now?: () => ms }

import {
    DEG, countAttachedEnds, attachmentAngles, matchAngleSets, checkEndpoints, maxTangentMisalignment,
    netCrossings, enclosedCharge, pointPolygonDistance, setCrossings, compareLinesInRect,
    compareContoursInRect, angleDiff, attachedEnds, attachmentAngle, pointPolylineDistance,
} from './fieldline-checks.js';
import { rng, rotateConfig, totalCharge } from './fieldline-configs.js';

const LINES_PER_UNIT = 4; // decision D4; env.models.LINES_PER_UNIT_CHARGE is checked in the first rule

export const viewOf = cfg => ({ width: cfg.width, height: cfg.height });

function model(env, kind, cfg) {
    return env.models.makeFieldModel(kind, cfg);
}
const electric = (env, cfg) => model(env, env.models.LINE_CHARGE_KIND, cfg);
const wire = (env, cfg) => model(env, env.models.WIRE_KIND, cfg);

const hasField = (cfg) => cfg.charges.some(c => !c.isTestCharge) || cfg.uniformX !== 0 || cfg.uniformY !== 0;
const fmt = (n, k = 1) => Number(n).toFixed(k);
const deg = r => r / DEG;

/** Trace the electric lines of a configuration; reports (and returns null for) a non-array result. */
function traceElectric(env, cfg, view, problems) {
    const m = electric(env, cfg);
    const lines = env.lines2d.electricFieldLines(m, view);
    if (!Array.isArray(lines)) {
        problems.push(`${cfg.name}: electricFieldLines did not return an array`);
        return { m, lines: null };
    }
    if (hasField(cfg) && lines.length === 0)
        problems.push(`${cfg.name}: electricFieldLines returned no lines for a model with ${m.sources.length} source(s) (is it still a stub?)`);
    return { m, lines: lines.length ? lines : null };
}

// ---------------------------------------------------------------------------

/** Each source of charge q has exactly 4|q| line ends attached (within 1 px). */
export function ruleExactCounts(env, configs) {
    const problems = [];
    if (env.models.LINES_PER_UNIT_CHARGE !== LINES_PER_UNIT)
        problems.push(`LINES_PER_UNIT_CHARGE is ${env.models.LINES_PER_UNIT_CHARGE}, tests assume ${LINES_PER_UNIT} (D4)`);
    for (const cfg of configs) {
        const { m, lines } = traceElectric(env, cfg, viewOf(cfg), problems);
        if (!lines) continue;
        const counts = countAttachedEnds(lines, m.sources, 1);
        m.sources.forEach((s, i) => {
            const expected = Math.round(LINES_PER_UNIT * Math.abs(s.charge));
            if (counts[i] !== expected)
                problems.push(`${cfg.name}: source ${i} (q=${s.charge} at ${fmt(s.x)},${fmt(s.y)}) has ${counts[i]} attached line ends, expected ${expected}`);
        });
        lines.forEach((l, li) => {
            if (typeof l.from === 'number' && l.from === l.to)
                problems.push(`${cfg.name}: line ${li} starts and ends on the same source ${l.from}`);
        });
    }
    return problems;
}

/** Every end is on a source or labelled 'infinity' and far away; from/to agree with the geometry. */
export function ruleEndpoints(env, configs) {
    const problems = [];
    for (const cfg of configs) {
        const view = viewOf(cfg);
        const { m, lines } = traceElectric(env, cfg, view, problems);
        if (!lines) continue;
        for (const p of checkEndpoints(lines, m.sources, view, { tol: 1, farFactor: 3 })) problems.push(`${cfg.name}: ${p}`);
        // oriented from + to -: a line must start on a positive (or at infinity) and end on a negative (or at infinity)
        lines.forEach((l, li) => {
            if (typeof l.from === 'number' && m.sources[l.from] && m.sources[l.from].charge <= 0)
                problems.push(`${cfg.name}: line ${li} starts on non-positive source ${l.from}`);
            if (typeof l.to === 'number' && m.sources[l.to] && m.sources[l.to].charge >= 0)
                problems.push(`${cfg.name}: line ${li} ends on non-negative source ${l.to}`);
        });
    }
    return problems;
}

/** Lines are tangent to model.fieldAt (oriented + -> -) at points outside the drawn charges (> 12 px from any source;
 * inside, the softened fieldAt and the exact stream function legitimately differ by a few degrees). */
export function ruleTangency(env, configs, { maxDeg = 3 } = {}) {
    const problems = [];
    for (const cfg of configs) {
        const view = viewOf(cfg);
        const { m, lines } = traceElectric(env, cfg, view, problems);
        if (!lines) continue;
        const cx = view.width / 2, cy = view.height / 2, diag = Math.hypot(view.width, view.height);
        let worst = { max: 0 }, worstLine = -1, checked = 0;
        lines.forEach((l, li) => {
            const r = maxTangentMisalignment(l.points, m.fieldAt, {
                sources: m.sources, minDist: 12, minField: m.K / 400, oriented: true,
                include: (x, y) => Math.hypot(x - cx, y - cy) < 2 * diag,
            });
            checked += r.n;
            if (r.max > worst.max) { worst = r; worstLine = li; }
        });
        if (deg(worst.max) >= maxDeg)
            problems.push(`${cfg.name}: line ${worstLine} (level ${lines[worstLine].level}) is ${fmt(deg(worst.max))} deg from the field direction (limit ${maxDeg}) near (${worst.at.map(v => fmt(v)).join(', ')})`);
        if (checked === 0) problems.push(`${cfg.name}: no segments were available to check tangency`);
    }
    return problems;
}

/** The same model traced with view x2 gives the same lines inside the original view (1 px; the view may set the
 * step size, which moves polyline vertices by a fraction of a pixel). */
export function ruleViewportIndependence(env, configs, { tol = 1 } = {}) {
    const problems = [];
    for (const cfg of configs) {
        const v1 = viewOf(cfg), v2 = { width: 2 * cfg.width, height: 2 * cfg.height };
        const p1 = [], p2 = [];
        const a = traceElectric(env, cfg, v1, p1), b = traceElectric(env, cfg, v2, p2);
        problems.push(...p1.map(p => `view 1x: ${p}`), ...p2.map(p => `view 2x: ${p}`));
        if (!a.lines || !b.lines) continue;
        // (a larger view legitimately contains more lines, e.g. of a uniform field, so only lines
        // passing through the original view are compared, in both directions)
        const rect = { x0: 0, y0: 0, x1: cfg.width, y1: cfg.height };
        for (const p of compareLinesInRect(a.lines, b.lines, rect, { tol })) problems.push(`${cfg.name}: ${p}`);
    }
    return problems;
}

/** Same for potential contours and magnetic lines ("view may only be used to choose resolution"). */
export function ruleViewportIndependenceContours(env, electricConfigs, magneticConfigs, { tol = 1 } = {}) {
    const problems = [];
    const run = (cfg, label, fn) => {
        const v1 = viewOf(cfg), v2 = { width: 2 * cfg.width, height: 2 * cfg.height };
        const a = fn(v1), b = fn(v2);
        if (!Array.isArray(a) || !Array.isArray(b)) { problems.push(`${cfg.name}: ${label} did not return an array`); return; }
        if (a.length === 0) { problems.push(`${cfg.name}: ${label} returned no contours`); return; }
        const rect = { x0: 0, y0: 0, x1: cfg.width, y1: cfg.height };
        for (const p of compareContoursInRect(a, b, rect, { tol })) problems.push(`${cfg.name} ${label}: ${p}`);
    };
    for (const cfg of electricConfigs) {
        const m = electric(env, cfg);
        if (hasField(cfg)) run(cfg, 'potentialContours', v => env.lines2d.potentialContours(m, v));
    }
    for (const cfg of magneticConfigs) {
        const m = wire(env, cfg);
        if (hasField(cfg)) run(cfg, 'magneticFieldLines', v => env.lines2d.magneticFieldLines(m, v));
    }
    return problems;
}

/**
 * Drag stability: move one charge in 1 px steps. Counts never change and attachment-angle sets at
 * every charge change by at most `maxStepDeg` per step (charges >= 30 px from all others).
 */
export function ruleDragStability(env, dragConfigs, { maxStepDeg = 3 } = {}) {
    const problems = [];
    for (const base of dragConfigs) {
        const { index, dx, dy } = base.move;
        const steps = Math.round(Math.max(Math.abs(dx), Math.abs(dy)));
        let prev = null, prevCounts = null;
        let worst = { step: -1, deg: 0, charge: -1 };
        let reported = 0;
        for (let s = 0; s <= steps; s++) {
            const charges = base.charges.map((c, i) => i === index ? { ...c, x: c.x + (dx * s) / steps, y: c.y + (dy * s) / steps } : { ...c });
            const moving = charges[index];
            const nearest = Math.min(...charges.filter((_, i) => i !== index).map(c => Math.hypot(c.x - moving.x, c.y - moving.y)));
            if (nearest < 30) { problems.push(`${base.name}: test setup error, moving charge comes within ${fmt(nearest)} px of another`); break; }
            const cfg = { ...base, charges };
            const p = [];
            const { m, lines } = traceElectric(env, cfg, viewOf(cfg), p);
            if (!lines) { if (s === 0) problems.push(...p); break; }
            const counts = countAttachedEnds(lines, m.sources, 1);
            const angles = attachmentAngles(lines, m.sources, { d: 15, tol: 1 });
            if (prev) {
                m.sources.forEach((src, i) => {
                    if (counts[i] !== prevCounts[i] && reported++ < 5)
                        problems.push(`${base.name}: step ${s}: source ${i} count changed from ${prevCounts[i]} to ${counts[i]}`);
                    const r = matchAngleSets(prev[i], angles[i]);
                    if (!r.ok) return;
                    if (deg(r.maxDiff) > worst.deg) worst = { step: s, deg: deg(r.maxDiff), charge: i };
                });
            }
            prev = angles; prevCounts = counts;
        }
        if (worst.deg > maxStepDeg)
            problems.push(`${base.name}: attachment angles at charge ${worst.charge} jump by ${fmt(worst.deg)} deg in one 1 px step (at step ${worst.step}; limit ${maxStepDeg})`);
    }
    return problems;
}

/** Neutral configurations rotated by 90 degrees about the centre have attachment angles rotated by 90 degrees. */
export function ruleRotationInvariance(env, configs, { tolDeg = 0.5 } = {}) {
    const problems = [];
    for (const cfg of configs) {
        if (Math.abs(totalCharge(cfg)) > 1e-9) { problems.push(`${cfg.name}: test setup error, configuration is not neutral`); continue; }
        const rotated = rotateConfig(cfg, 1);
        const a = traceElectric(env, cfg, viewOf(cfg), problems), b = traceElectric(env, rotated, viewOf(rotated), problems);
        if (!a.lines || !b.lines) continue;
        const A = attachmentAngles(a.lines, a.m.sources), B = attachmentAngles(b.lines, b.m.sources);
        a.m.sources.forEach((s, i) => {
            const r = matchAngleSets(A[i], B[i], Math.PI / 2);
            if (!r.ok) problems.push(`${cfg.name}: source ${i} has ${r.sizeA} attached lines but ${r.sizeB} after rotating`);
            else if (deg(r.maxDiff) > tolDeg)
                problems.push(`${cfg.name}: attachment angles at source ${i} differ by ${fmt(deg(r.maxDiff), 2)} deg from the rotated originals (limit ${tolDeg}); before ${A[i].map(x => fmt(deg(x))).join(' ')}, after ${B[i].map(x => fmt(deg(x))).join(' ')}`);
        });
    }
    return problems;
}

/** A +2/-2 horizontal pair has a straight line along the axis; a +2/+2 pair has outward axis lines. */
export function ruleAxisLine(env, axis, { tolDeg = 0.5 } = {}) {
    const problems = [];
    const check = (cfg, plusIndex, wantAngle, what, straightTo = null) => {
        const { m, lines } = traceElectric(env, cfg, viewOf(cfg), problems);
        if (!lines) return;
        const S = m.sources[plusIndex];
        const ends = attachedEnds(lines, m.sources, 1)[plusIndex];
        const found = ends.map(e => ({ ...e, ang: attachmentAngle(lines[e.line].points, e.end, S, 15) }))
            .filter(e => Math.abs(angleDiff(e.ang, wantAngle)) <= tolDeg * DEG);
        if (found.length === 0) {
            problems.push(`${cfg.name}: no line leaves source ${plusIndex} ${what} (wanted attachment angle ${fmt(deg(wantAngle))} deg +- ${tolDeg}); angles are ${ends.map(e => fmt(deg(attachmentAngle(lines[e.line].points, e.end, S, 15)))).join(' ')}`);
            return;
        }
        if (straightTo) {
            const axisLine = [[S.x, S.y], straightTo];
            const l = lines[found[0].line];
            const dev = Math.max(...l.points.map(p => pointPolylineDistance(p, axisLine)));
            if (dev > 0.5) problems.push(`${cfg.name}: the axis line ${what} deviates from the axis by ${fmt(dev, 2)} px (limit 0.5)`);
        }
    };
    const u = axis.unlike, l = axis.like;
    // unlike: + at index 0, - at index 1 on the right: towards the minus charge is angle 0 (+x)
    check(u, 0, 0, 'towards the negative charge', [u.charges[1].x, u.charges[1].y]);
    // like: charge 0 is on the left, so outwards is +pi; charge 1 outwards is 0
    check(l, 0, Math.PI, 'outwards along the axis');
    check(l, 1, 0, 'outwards along the axis');
    return problems;
}

/** Potential contours cross electric field lines at 90 +- tolDeg degrees. */
export function rulePerpendicularContours(env, configs, { tolDeg = 6, minCrossings = 4 } = {}) {
    // tolerance allows for measuring with polyline chords on tightly curved contours near charges
    const problems = [];
    for (const cfg of configs) {
        const view = viewOf(cfg);
        const { m, lines } = traceElectric(env, cfg, view, problems);
        if (!lines) continue;
        const contours = env.lines2d.potentialContours(m, view);
        if (!Array.isArray(contours) || contours.length === 0) {
            problems.push(`${cfg.name}: potentialContours returned no contours`);
            continue;
        }
        const rect = { x0: 0, y0: 0, x1: cfg.width, y1: cfg.height };
        let n = 0, worst = null;
        for (const c of setCrossings(contours, lines, rect)) {
            const f = m.fieldAt(c.x, c.y);
            if (Math.hypot(f.u, f.v) < m.K / 400) continue; // next to a saddle point
            if (m.sources.some(s => Math.hypot(c.x - s.x, c.y - s.y) < 8)) continue;
            n++;
            const off = 90 - deg(c.angle);
            if (!worst || off > worst.off) worst = { off, x: c.x, y: c.y, level: contours[c.a].level };
        }
        if (m.sources.length && n < minCrossings) problems.push(`${cfg.name}: only ${n} contour/field-line crossings inside the view`);
        else if (worst && worst.off > tolDeg)
            problems.push(`${cfg.name}: a potential contour (level ${fmt(worst.level, 3)}) crosses a field line at ${fmt(90 - worst.off)} deg (limit 90 +- ${tolDeg}) near (${fmt(worst.x)}, ${fmt(worst.y)})`);
    }
    return problems;
}

/** Random closed polygons, avoiding sources by >= 15 px, have net signed crossings = 4 x enclosed charge. */
export function ruleGauss(env, configs, { seed = 99, polygonsPerConfig = 8, clearance = 15 } = {}) {
    const problems = [];
    configs.forEach((cfg, ci) => {
        const { m, lines } = traceElectric(env, cfg, viewOf(cfg), problems);
        if (!lines) return;
        const rand = rng(seed * 1000 + ci);
        let made = 0, attempts = 0, enclosingSome = 0;
        while (made < polygonsPerConfig && attempts++ < 2000) {
            // alternate between polygons centred near a source and anywhere
            let cx, cy;
            if (m.sources.length && made % 2 === 0) {
                const s = m.sources[Math.floor(rand() * m.sources.length)];
                cx = s.x + (rand() - 0.5) * 60; cy = s.y + (rand() - 0.5) * 60;
            } else { cx = rand() * cfg.width; cy = rand() * cfg.height; }
            const n = 5 + Math.floor(rand() * 6), base = 25 + rand() * 150;
            const poly = [];
            for (let i = 0; i < n; i++) {
                const t = ((i + 0.6 * (rand() - 0.5)) / n) * 2 * Math.PI, r = base * (0.55 + 0.7 * rand());
                poly.push([cx + r * Math.cos(t), cy + r * Math.sin(t)]);
            }
            if (m.sources.some(s => pointPolygonDistance([s.x, s.y], poly) < clearance)) continue;
            made++;
            const q = enclosedCharge(m.sources, poly);
            if (q !== 0) enclosingSome++;
            const got = netCrossings(lines, poly);
            const expected = LINES_PER_UNIT * q;
            if (got.net !== expected)
                problems.push(`${cfg.name}: polygon ${JSON.stringify(poly.map(p => p.map(v => Math.round(v))))} encloses charge ${q}: net crossings ${got.net} (out ${got.out}, in ${got.in}), expected ${expected}`);
        }
        if (made < polygonsPerConfig) problems.push(`${cfg.name}: test setup error, only generated ${made} polygons`);
    });
    return problems;
}

/** Magnetic contours: closed without a uniform field, tangent to B, levels spaced by contourStep. */
export function ruleMagnetic(env, configs, { tangentDeg = 3 } = {}) {
    const problems = [];
    for (const cfg of configs) {
        const view = viewOf(cfg);
        const m = wire(env, cfg);
        if (!hasField(cfg)) continue;
        const lines = env.lines2d.magneticFieldLines(m, view);
        if (!Array.isArray(lines) || lines.length === 0) {
            problems.push(`${cfg.name}: magneticFieldLines returned no contours for ${m.sources.length} source(s) (is it still a stub?)`);
            continue;
        }
        const step = env.lines2d.contourStep(m);
        if (!(step > 0) || !isFinite(step)) { problems.push(`${cfg.name}: contourStep is ${step}, expected a positive number`); continue; }
        // radius of the drawn wire disk, inside which contours may be omitted
        const diskRadius = env.lines2d.SOURCE_DISK_RADIUS ?? 10;

        // closed loops when there is no uniform field. A neutral set of wires also has an
        // infinite zero-level contour (and arbitrarily large loops at tiny levels): those are
        // allowed to be open, but only if they extend far outside the view
        const noUniform = cfg.uniformX === 0 && cfg.uniformY === 0;
        const neutral = Math.abs(totalCharge(cfg)) < 1e-9;
        if (noUniform) {
            lines.forEach((l, li) => {
                if (l.closed === true) return;
                // clipped where it enters the drawn disk of a wire (both ends on a disk)
                const onDisk = p => m.sources.some(src => Math.abs(Math.hypot(p[0] - src.x, p[1] - src.y) - diskRadius) <= 1);
                const clipped = onDisk(l.points[0]) && onDisk(l.points[l.points.length - 1]);
                const exempt = neutral || clipped;
                if (!exempt) problems.push(`${cfg.name}: contour ${li} (level ${fmt(l.level, 3)}) is open although there is no uniform field and the net current is not zero`);
            });
        }
        for (const p of checkEndpoints(lines, m.sources, view, { labelled: false, diskRadius })) problems.push(`${cfg.name}: ${p}`);

        // tangency to B (undirected)
        const cx = view.width / 2, cy = view.height / 2, diag = Math.hypot(view.width, view.height);
        let worst = { max: 0 }, worstLine = -1;
        lines.forEach((l, li) => {
            const r = maxTangentMisalignment(l.points, m.fieldAt, {
                sources: m.sources, minDist: 12, minField: m.K / 400, oriented: false,
                include: (x, y) => Math.hypot(x - cx, y - cy) < 2 * diag,
            });
            if (r.max > worst.max) { worst = r; worstLine = li; }
        });
        if (deg(worst.max) >= tangentDeg)
            problems.push(`${cfg.name}: contour ${worstLine} (level ${fmt(lines[worstLine].level, 3)}) is ${fmt(deg(worst.max))} deg from the field direction near (${worst.at.map(v => fmt(v)).join(', ')}) (limit ${tangentDeg})`);

        // level spacing: every level is an integer multiple of contourStep ...
        const levels = [...new Set(lines.map(l => l.level))];
        for (const lv of levels) {
            const k = lv / step;
            if (Math.abs(k - Math.round(k)) > 1e-6 * Math.max(1, Math.abs(k)))
                problems.push(`${cfg.name}: contour level ${lv} is not an integer multiple of contourStep ${step}`);
        }
        // ... and every multiple within the range of A_z over the view is present (ignoring the
        // 15 px around each wire, where A_z diverges logarithmically and tiny loops may be dropped)
        let lo = Infinity, hi = -Infinity;
        for (let x = 0; x <= cfg.width; x += 4) for (let y = 0; y <= cfg.height; y += 4) {
            if (m.sources.some(s => Math.hypot(x - s.x, y - s.y) < 15)) continue;
            const a = m.potentialAt(x, y);
            lo = Math.min(lo, a); hi = Math.max(hi, a);
        }
        const have = new Set(levels.map(lv => Math.round(lv / step)));
        const missing = [];
        for (let k = Math.ceil(lo / step + 0.05); k <= Math.floor(hi / step - 0.05); k++)
            if (k !== 0 && !have.has(k)) missing.push(k);
        if (missing.length)
            problems.push(`${cfg.name}: no contour for level(s) ${missing.slice(0, 6).map(k => k + ' x step').join(', ')}${missing.length > 6 ? ' ...' : ''} although A_z takes those values inside the view (range ${fmt(lo / step, 2)}..${fmt(hi / step, 2)} steps)`);
    }
    return problems;
}

/** contourStep is proportional to model.K: the same relative picture on a different canvas size. */
export function ruleContourStepScaling(env, configs) {
    const problems = [];
    for (const cfg of configs) {
        const mA = electric(env, cfg), s = 2;
        const big = { ...cfg, width: cfg.width * s, height: cfg.height * s, charges: cfg.charges.map(c => ({ ...c, x: c.x * s, y: c.y * s })) };
        const mB = electric(env, big);
        const a = env.lines2d.contourStep(mA) / mA.K, b = env.lines2d.contourStep(mB) / mB.K;
        if (!(a > 0)) problems.push(`${cfg.name}: contourStep must be positive`);
        else if (Math.abs(a - b) > 1e-9 * a) problems.push(`${cfg.name}: contourStep / K is ${a} at ${cfg.width} px but ${b} at ${big.width} px (the picture must not depend on canvas size)`);
    }
    return problems;
}

/** Median time of 10 calls (ms) under the budget for 4-charge configurations. */
export function ruleTimeBudget(env, configs, { budgetMs = 30, calls = 10 } = {}) {
    const problems = [];
    const now = env.now || (() => performance.now());
    const median = a => [...a].sort((x, y) => x - y)[Math.floor(a.length / 2)];
    for (const cfg of configs) {
        const view = viewOf(cfg);
        const me = electric(env, cfg), mw = wire(env, cfg);
        const tE = [], tM = [];
        let nE = 0, nM = 0;
        for (let i = 0; i < calls; i++) {
            let t0 = now();
            nE = env.lines2d.electricFieldLines(me, view).length + env.lines2d.potentialContours(me, view).length;
            tE.push(now() - t0);
            t0 = now();
            nM = env.lines2d.magneticFieldLines(mw, view).length;
            tM.push(now() - t0);
        }
        if (nE === 0) problems.push(`${cfg.name}: electric lines and potential contours are empty (is it still a stub?)`);
        if (nM === 0) problems.push(`${cfg.name}: magnetic lines are empty (is it still a stub?)`);
        if (median(tE) >= budgetMs) problems.push(`${cfg.name}: electricFieldLines + potentialContours median ${fmt(median(tE), 1)} ms (budget ${budgetMs} ms)`);
        if (median(tM) >= budgetMs) problems.push(`${cfg.name}: magneticFieldLines median ${fmt(median(tM), 1)} ms (budget ${budgetMs} ms)`);
    }
    return problems;
}
