// Rule tests for field lines of the grid solvers (electrostatic_fourier, dynamic): gridfield.js +
// fieldlines2d.js on top of the WASM stream function (maxwell/src/streamfn.rs). Issue #3, stage 2b.
//
// The grid world ends at the edge of the canvas, so (unlike the analytic models) lines end on a
// charge, on the edge of the canvas, or where they fade out in a region with (numerically) no
// field at all, e.g. where radiation has not arrived yet and the static fields cancel.
//
// Scenarios: the lecture examples that use the grid solvers, random 3-4 charge configurations
// (Fourier), a net charge (Fourier), and the dynamic solver after dragging a charge and after
// shaking one (radiation, with closed field lines).

import { test, expect } from '@playwright/test';
import { examples } from './helpers.js';
import { lectureConfigs, randomConfigs } from './fieldline-configs.js';

const fourier = [
    ...lectureConfigs(examples, 'electrostatic_fourier'),
    ...randomConfigs(6, { seed: 777 }),
    { name: 'net-charge', charges: [{ x: 300, y: 150, charge: 2 }], uniformX: 0, uniformY: 0, width: 668, height: 335 },
].map(c => ({ ...c, solver: 'electrostatic_fourier' }));

const l14 = lectureConfigs(examples, 'dynamic').find(c => c.name === 'lecture14-2');
const l16 = lectureConfigs(examples, 'dynamic').find(c => c.name === 'lecture16-1');
const dynamic = [
    { ...l14, name: 'lecture14-2 static', solver: 'dynamic' },
    // drag the +3 charge 60 px right and 20 px down over 0.3 s, then wait 0.5 s
    { ...l14, name: 'lecture14-2 dragged', solver: 'dynamic', motion: { index: 0, kind: 'drag', dx: 60, dy: 20, ms: 300, waitMs: 500 } },
    // shake the +2 charge vertically (+-25 px, period 0.5 s) for 1 s, then wait 0.4 s
    { ...l16, name: 'lecture16-1 shaken', solver: 'dynamic', motion: { index: 0, kind: 'shake', amp: 25, periodMs: 500, ms: 1000, waitMs: 400 } },
];

/** Compute the grid field lines for each config in the page; return per-config data for the checks. */
async function gridLines(page, configs) {
    await page.goto('/tests/blank.html');
    return page.evaluate(async (configs) => {
        const wasm = await import('/maxwell/out/maxwell.js');
        await wasm.default();
        const { makeGridFieldModel } = await import('/gridfield.js');
        const { electricFieldLines } = await import('/fieldlines2d.js');
        const out = [];
        for (const cfg of configs) {
            const { width: W, height: H } = cfg;
            const field = new wasm.FieldConfiguration(W, H, 111, 111);
            let charges = cfg.charges.map(c => ({ ...c }));
            const sources = () => charges.filter(c => !c.isTestCharge);
            field.set_charges(sources());
            field.set_uniform_field(cfg.uniformX, cfg.uniformY);
            if (cfg.solver === 'dynamic') {
                // as ui.js: ticks of 0.1 cell, internal time = ms / 5, 60 frames per second
                const step = 0.1 * W / 111, frameMs = 1000 / 60;
                const m = cfg.motion;
                const start = { ...charges[m?.index ?? 0] };
                const frames = m ? Math.round((m.ms + m.waitMs) / frameMs) : 0;
                for (let f = 1; f <= frames; f++) {
                    const t = f * frameMs;
                    if (t <= m.ms) {
                        const c = charges[m.index];
                        if (m.kind === 'drag') { c.x = start.x + m.dx * t / m.ms; c.y = start.y + m.dy * t / m.ms; }
                        else { c.y = start.y + m.amp * Math.sin(2 * Math.PI * t / m.periodMs); }
                        if (t + frameMs > m.ms) { c.x = m.kind === 'drag' ? start.x + m.dx : start.x; c.y = m.kind === 'drag' ? start.y + m.dy : start.y; }
                        field.set_charges(sources());
                    }
                    let dt = frameMs / 5;
                    while (dt > 0) { const s = Math.min(step, dt); field.tick(s); dt -= s; }
                }
            } else {
                field.reset_fields();
            }
            const model = makeGridFieldModel({ field, charges, uniformX: cfg.uniformX, uniformY: cfg.uniformY, width: W, height: H });
            const lines = electricFieldLines(model, { width: W, height: H });
            // the field (from psi) at every line end and at segment midpoints, for the checks below
            const fieldAtEnds = lines.map(l => [l.points[0], l.points[l.points.length - 1]].map(p => model.fieldAt(p[0], p[1])));
            const times = [];
            for (let k = 0; k < 5; k++) {
                const t0 = performance.now();
                electricFieldLines(makeGridFieldModel({ field, charges, uniformX: cfg.uniformX, uniformY: cfg.uniformY, width: W, height: H }), { width: W, height: H });
                times.push(performance.now() - t0);
            }
            times.sort((a, b) => a - b);
            const tangency = [];
            for (const l of lines) {
                let worst = 0, at = null;
                for (let i = 0; i + 1 < l.points.length; i++) {
                    const [p, q] = [l.points[i], l.points[i + 1]];
                    const mx = (p[0] + q[0]) / 2, my = (p[1] + q[1]) / 2;
                    if (model.sources.some(s => Math.hypot(mx - s.x, my - s.y) < 12)) continue;
                    // a final segment cut off at the edge of the canvas is an off-centre piece of a chord
                    const edge = (r) => Math.min(r[0], W - r[0], r[1], H - r[1]) < 1e-6;
                    if (edge(p) || edge(q)) continue;
                    const f = model.fieldAt(mx, my), dx = q[0] - p[0], dy = q[1] - p[1];
                    if (Math.hypot(f.u, f.v) < 0.05 || Math.hypot(dx, dy) === 0) continue;
                    const mis = Math.atan2(Math.abs(dx * f.v - dy * f.u), dx * f.u + dy * f.v);
                    if (mis > worst) { worst = mis; at = [mx, my]; }
                }
                tangency.push({ worst, at });
            }
            out.push({
                name: cfg.name, width: W, height: H,
                sources: model.sources.map(s => ({ x: s.x, y: s.y, charge: s.charge })),
                lines: lines.map(l => ({ points: l.points, from: l.from, to: l.to, closed: !!l.closed })),
                fieldAtEnds, tangency, medianMs: times[2],
            });
        }
        return out;
    }, configs);
}

const onEdge = (p, W, H) => Math.min(Math.abs(p[0]), Math.abs(p[0] - W), Math.abs(p[1]), Math.abs(p[1] - H)) < 1e-6;
const onSource = (p, sources) => sources.some(s => Math.hypot(p[0] - s.x, p[1] - s.y) <= 1);
const fmt = (v, k = 1) => Number(v).toFixed(k);

function checks(results, { tangentDeg = 3 } = {}) {
    const problems = [];
    for (const r of results) {
        // 1. exactly 4|q| line ends on each source
        r.sources.forEach((s, i) => {
            let n = 0;
            for (const l of r.lines) for (const p of [l.points[0], l.points[l.points.length - 1]])
                if (!l.closed && Math.hypot(p[0] - s.x, p[1] - s.y) <= 1) n++;
            if (n !== Math.round(4 * Math.abs(s.charge)))
                problems.push(`${r.name}: source ${i} (q = ${s.charge}) has ${n} line ends, expected ${4 * Math.abs(s.charge)}`);
        });
        r.lines.forEach((l, li) => {
            const ends = [l.points[0], l.points[l.points.length - 1]];
            if (l.closed) {
                // 2a. closed lines close up
                const gap = Math.hypot(ends[0][0] - ends[1][0], ends[0][1] - ends[1][1]);
                if (gap > 20) problems.push(`${r.name}: closed line ${li} has a gap of ${fmt(gap)} px`);
                return;
            }
            // 2b. open lines end on a source, on the edge of the canvas, or fade out where there is no field
            ends.forEach((p, k) => {
                const f = r.fieldAtEnds[li][k];
                if (!onSource(p, r.sources) && !onEdge(p, r.width, r.height) && Math.hypot(f.u, f.v) > 0.05)
                    problems.push(`${r.name}: line ${li} ends at (${fmt(p[0])}, ${fmt(p[1])}), not on a source or the edge, where |E| = ${fmt(Math.hypot(f.u, f.v), 3)}`);
            });
            // 3. tangent to the field (the field from psi; away from the sources)
            const t = r.tangency[li];
            if (t.worst > tangentDeg * Math.PI / 180)
                problems.push(`${r.name}: line ${li} is ${fmt(t.worst * 180 / Math.PI)} deg off the field direction at (${fmt(t.at[0])}, ${fmt(t.at[1])})`);
        });
    }
    return problems;
}

function gauss(results, { seed = 5, polygons = 8, clearance = 15 } = {}) {
    let s = seed;
    const rand = () => (s = (s * 16807) % 2147483647) / 2147483647;
    const crossings = (points, poly) => {
        // signed crossings, outward positive (polygon orientation from its signed area)
        let a2 = 0;
        for (let i = 0; i < poly.length; i++) { const p = poly[i], q = poly[(i + 1) % poly.length]; a2 += p[0] * q[1] - q[0] * p[1]; }
        const orient = a2 > 0 ? 1 : -1;
        let net = 0;
        for (let i = 0; i + 1 < points.length; i++) {
            const p = points[i], q = points[i + 1], rx = q[0] - p[0], ry = q[1] - p[1];
            for (let j = 0; j < poly.length; j++) {
                const a = poly[j], b = poly[(j + 1) % poly.length], sx = b[0] - a[0], sy = b[1] - a[1];
                const den = rx * sy - ry * sx;
                if (den === 0) continue;
                const t = ((a[0] - p[0]) * sy - (a[1] - p[1]) * sx) / den, u = ((a[0] - p[0]) * ry - (a[1] - p[1]) * rx) / den;
                if (t < 0 || t >= 1 || u < 0 || u >= 1) continue;
                net += orient * den > 0 ? 1 : -1;
            }
        }
        return net;
    };
    const inside = (p, poly) => {
        let c = false;
        for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
            const [xi, yi] = poly[i], [xj, yj] = poly[j];
            if ((yi > p[1]) !== (yj > p[1]) && p[0] < (xj - xi) * (p[1] - yi) / (yj - yi) + xi) c = !c;
        }
        return c;
    };
    const segDist = (p, a, b) => {
        const dx = b[0] - a[0], dy = b[1] - a[1], l2 = dx * dx + dy * dy;
        const t = l2 ? Math.max(0, Math.min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2)) : 0;
        return Math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy);
    };
    const problems = [];
    for (const r of results) {
        let made = 0, attempts = 0;
        while (made < polygons && attempts++ < 2000) {
            const src = r.sources.length && made % 2 === 0 ? r.sources[Math.floor(rand() * r.sources.length)] : null;
            const cx = src ? src.x + (rand() - 0.5) * 60 : rand() * r.width, cy = src ? src.y + (rand() - 0.5) * 60 : rand() * r.height;
            const n = 5 + Math.floor(rand() * 6), base = 25 + rand() * 120, poly = [];
            for (let i = 0; i < n; i++) {
                const t = ((i + 0.6 * (rand() - 0.5)) / n) * 2 * Math.PI, rr = base * (0.55 + 0.7 * rand());
                poly.push([cx + rr * Math.cos(t), cy + rr * Math.sin(t)]);
            }
            // keep inside the canvas (the grid world) and clear of the sources
            if (poly.some(p => p[0] < 1 || p[1] < 1 || p[0] > r.width - 1 || p[1] > r.height - 1)) continue;
            if (r.sources.some(s => poly.some((p, i) => segDist([s.x, s.y], p, poly[(i + 1) % n]) < clearance))) continue;
            made++;
            const q = r.sources.filter(s => inside([s.x, s.y], poly)).reduce((t, s) => t + s.charge, 0);
            const net = r.lines.reduce((t, l) => t + crossings(l.closed ? [...l.points, l.points[0]] : l.points, poly), 0);
            if (net !== 4 * q) problems.push(`${r.name}: polygon around (${fmt(cx, 0)}, ${fmt(cy, 0)}) encloses charge ${q}: net crossings ${net}, expected ${4 * q}`);
        }
    }
    return problems;
}

function report(problems, title) {
    if (problems.length) throw new Error(`${title}: ${problems.length} violation(s)\n${problems.slice(0, 12).map(p => `  - ${p}`).join('\n')}`);
    expect(problems).toEqual([]);
}

// The dynamic field has cell-scale structure (e.g. next to the absorbing layer the direction can
// swing by ~40 deg over two cells), so a chord within the 0.1 px sagitta tolerance can be a few
// degrees off the local direction.
const TANGENT_DEG = { 'Fourier solver': 3, 'dynamic solver': 8 };

for (const [group, configs] of [['Fourier solver', fourier], ['dynamic solver', dynamic]]) {
    test.describe(`grid field lines: ${group}`, () => {
        let results;
        test.beforeEach(async ({ page }) => { results ??= await gridLines(page, configs); });

        test(`4|q| line ends per charge; lines end on a charge, the edge, or where there is no field; tangent to the field (< ${TANGENT_DEG[group]} deg)`, () => {
            report(checks(results, { tangentDeg: TANGENT_DEG[group] }), 'Grid field lines');
        });
        test("Gauss's law: net crossings of random polygons = 4 x enclosed charge", () => {
            report(gauss(results), "Gauss's law");
        });
        test('time budget: median < 30 ms per configuration', () => {
            report(results.filter(r => r.medianMs > 30).map(r => `${r.name}: ${fmt(r.medianMs)} ms`), 'Time budget');
        });
    });
}

test('the shaken dipole radiates closed field lines', async ({ page }) => {
    const [r] = await gridLines(page, [dynamic[2]]);
    expect(r.lines.filter(l => l.closed).length).toBeGreaterThan(0);
});

// Causality (dynamic solver): while a charge is dragged, the field lines ahead of the wavefront must
// not move: by at most 0.5 px, or 2% of the local line spacing where the field is weak. Light travels
// 1 px per 5 ms; the margin covers the smeared charge (8 cells, ~63 px in x).
test('dynamic solver: field lines outside the light cone do not move during a drag', async ({ page }) => {
    await page.goto('/tests/blank.html');
    const problems = await page.evaluate(async () => {
        const wasm = await import('/maxwell/out/maxwell.js');
        await wasm.default();
        const { makeGridFieldModel } = await import('/gridfield.js');
        const { electricFieldLines } = await import('/fieldlines2d.js');
        const W = 668, H = 335, frameMs = 1000 / 60, step = 0.1 * W / 111;
        const out = [];
        for (const scenario of [
            { name: '+-3 pair, +3 dragged', charges: [{ x: 467.6, y: 167.5, charge: 3 }, { x: 200.4, y: 167.5, charge: -3 }], move: [60, 20] },
            { name: 'single +2 dragged', charges: [{ x: 250, y: 150, charge: 2 }], move: [-40, 30] },
        ]) {
            const field = new wasm.FieldConfiguration(W, H, 111, 111);
            const charges = scenario.charges.map(c => ({ ...c }));
            const start = { ...charges[0] };
            field.set_charges(charges);
            let model;
            const lines = () => {
                model = makeGridFieldModel({ field, charges, width: W, height: H });
                return electricFieldLines(model, { width: W, height: H });
            };
            const before = lines();
            // the local spacing of the field lines: psi changes by 1/4 between neighbours
            const spacing = (x, y) => { const e = model.fieldAt(x, y); return 2 * Math.PI * model.K / (4 * Math.hypot(e.u, e.v)); };
            // points of `a` farther than `r` from the moving charge's path, and their distance to the nearest segment of `b`
            const far = (a, b, r, u) => {
                const [ex, ey] = [scenario.move[0] * u, scenario.move[1] * u];
                const fromPath = (x, y) => {
                    const l2 = ex * ex + ey * ey, t = l2 ? Math.max(0, Math.min(1, ((x - start.x) * ex + (y - start.y) * ey) / l2)) : 0;
                    return Math.hypot(x - start.x - t * ex, y - start.y - t * ey);
                };
                // the worst ratio of the displacement to the allowed one: max(0.5 px, 2% of the line spacing)
                let worst = 0, at = null, n = 0, px = 0;
                for (const l of a) for (const [x, y] of l.points) {
                    if (fromPath(x, y) < r || x < 1 || y < 1 || x > W - 1 || y > H - 1) continue;
                    n++;
                    let best = Infinity;
                    for (const m of b) for (let k = 0; k + 1 < m.points.length; k++) {
                        const [p, q] = [m.points[k], m.points[k + 1]];
                        const dx = q[0] - p[0], dy = q[1] - p[1], l2 = dx * dx + dy * dy;
                        const t = l2 ? Math.max(0, Math.min(1, ((x - p[0]) * dx + (y - p[1]) * dy) / l2)) : 0;
                        best = Math.min(best, Math.hypot(x - p[0] - t * dx, y - p[1] - t * dy));
                    }
                    const ratio = best / Math.max(0.5, 0.02 * spacing(x, y));
                    if (ratio > worst) { worst = ratio; at = [x, y]; px = best; }
                }
                return { worst, at, n, px };
            };
            for (let f = 1; f <= 30; f++) {
                const t = f * frameMs, u = Math.min(1, t / 300);
                charges[0].x = start.x + scenario.move[0] * u;
                charges[0].y = start.y + scenario.move[1] * u;
                field.set_charges(charges);
                let dt = frameMs / 5;
                while (dt > 0) { const s = Math.min(step, dt); field.tick(s); dt -= s; }
                const now = lines();
                if (f % 5) continue;
                const r = t / 5 + 80;
                for (const [a, b, what] of [[now, before, 'a line point'], [before, now, 'an earlier line point']]) {
                    const { worst, at, n, px } = far(a, b, r, u);
                    if (n > 0 && worst > 1)
                        out.push(`${scenario.name}, ${t.toFixed(0)} ms: ${what} at (${at[0].toFixed(0)}, ${at[1].toFixed(0)}), outside the light cone (r > ${r.toFixed(0)} px), is ${px.toFixed(2)} px (${(100 * px / spacing(at[0], at[1])).toFixed(1)}% of the line spacing) from the other frame's lines`);
                }
            }
        }
        return out;
    });
    expect(problems).toEqual([]);
});

// As for the analytic models (decision D5): a symmetric +-3 pair has a field line along its axis,
// for the Fourier solver and for the dynamic solver's starting state.
test('grid solvers: a +-3 pair has a field line along its axis', async ({ page }) => {
    await page.goto('/tests/blank.html');
    const found = await page.evaluate(async () => {
        const wasm = await import('/maxwell/out/maxwell.js');
        await wasm.default();
        const { makeGridFieldModel } = await import('/gridfield.js');
        const { electricFieldLines } = await import('/fieldlines2d.js');
        const charges = [{ x: 467.6, y: 167.5, charge: 3 }, { x: 200.4, y: 167.5, charge: -3 }];
        const out = [];
        for (const dynamic of [false, true]) {
            const field = new wasm.FieldConfiguration(668, 335, 111, 111);
            field.set_charges(charges);
            if (!dynamic) field.reset_fields();
            const lines = electricFieldLines(makeGridFieldModel({ field, charges, width: 668, height: 335 }), { width: 668, height: 335 });
            out.push(lines.some(l => l.from === 0 && l.to === 1 && l.points.every(([, y]) => Math.abs(y - 167.5) < 0.5)));
        }
        return out;
    });
    expect(found).toEqual([true, true]);
});

// The same through the real UI (ui.js's frame loop): drag the +3 charge of lecture 14's dynamic
// example with field lines shown; in the left third of the canvas (the -3 charge, at least 200 px
// beyond the light front while the drag lasts) lines may only shift by a fraction of a pixel: no
// pixel's RGBA may change by more than 128 in total (a line moving by a whole pixel changes pixels
// by ~255). This caught the field being ticked with the previous frame's positions, which made the
// far field lines jitter by a pixel or more while a charge was dragged.
test('dynamic solver in the UI: far field lines stay still while a charge is dragged', async ({ page }) => {
    const { examples } = await import('./helpers.js');
    const ex = examples.find(e => e.id === 'lecture14-2');
    const state = { ...ex.state, plotType: 'fieldline' };
    await page.clock.install({ time: new Date('2025-01-01T00:00:00Z') });
    await page.clock.pauseAt(new Date('2025-01-01T00:00:01Z'));
    await page.goto('/tests/harness.html?attrs=' + encodeURIComponent(JSON.stringify({ meme: JSON.stringify(state) })));
    for (let i = 0; i < 40; i++) { await page.clock.runFor(48); await page.waitForTimeout(10); }
    const canvas = page.locator('meme-embed canvas').first();
    const leftThird = () => canvas.evaluate(c => Array.from(c.getContext('2d').getImageData(0, 0, Math.floor(c.width / 3), c.height).data));
    const before = await leftThird();
    const box = await canvas.boundingBox();
    const a = { x: box.x + 0.7 * box.width, y: box.y + 0.5 * box.height };
    await page.mouse.move(Math.round(a.x), Math.round(a.y));
    await page.mouse.down();
    const jumps = [];
    for (let i = 1; i <= 18; i++) {   // 300 ms, 60 x 20 canvas px, uneven steps as from a real mouse
        const u = (i + 0.4 * Math.sin(i)) / 18.4;
        await page.mouse.move(Math.round(a.x + 60 * u * box.width / 668), Math.round(a.y + 20 * u * box.height / 335));
        await page.clock.runFor(16.7);
        const now = await leftThird();
        let worst = 0;
        for (let k = 0; k < now.length; k += 4) {
            const d = Math.abs(now[k] - before[k]) + Math.abs(now[k + 1] - before[k + 1]) + Math.abs(now[k + 2] - before[k + 2]) + Math.abs(now[k + 3] - before[k + 3]);
            if (d > worst) worst = d;
        }
        if (worst > 128) jumps.push(`frame ${i}: a pixel changed by ${worst}`);
    }
    await page.mouse.up();
    expect(jumps).toEqual([]);
});
