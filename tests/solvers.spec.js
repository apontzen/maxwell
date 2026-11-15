// Unit tests of the solver registry (solvers.js) and the 2D field models (fieldmodels.js):
// in field-line mode every overlay must come from the same 2D model (issue #3).

import { test, expect } from './helpers.js';

test.beforeEach(async ({ page }) => {
    await page.goto('/tests/blank.html');
    // solvers.js imports the WASM module, which must be initialised before use
    await page.evaluate(async () => { const wasm = await import('/maxwell/out/maxwell.js'); await wasm.default(); });
});

test('effective solver and aliases', async ({ page }) => {
    const r = await page.evaluate(async () => {
        const s = await import('/solvers.js');
        return {
            alias: s.canonicalSolver('magnetostatic_direct'),
            e3fl: s.effectiveSolver('electrostatic_direct', 'fieldline'),
            e3q: s.effectiveSolver('electrostatic_direct', 'quiver'),
            m3fl: s.effectiveSolver('magnetostatic_direct', 'fieldline'),
            m3q: s.effectiveSolver('magnetostatic_direct', 'quiver'),
            e2: s.effectiveSolver('electrostatic_line_charge', 'quiver'),
            w2: s.effectiveSolver('magnetostatic_wire', 'quiver'),
            fourier: s.effectiveSolver('electrostatic_fourier', 'fieldline'),
            magnetic: ['magnetostatic_direct', 'magnetostatic_wire', 'electrostatic_direct'].map(s.isMagnetic),
        };
    });
    expect(r).toEqual({
        alias: 'magnetostatic_moving_charge',
        e3fl: 'electrostatic_line_charge', e3q: 'electrostatic_direct',
        m3fl: 'magnetostatic_wire', m3q: 'magnetostatic_moving_charge',
        e2: 'electrostatic_line_charge', w2: 'magnetostatic_wire',
        fourier: 'electrostatic_fourier',
        magnetic: [true, true, false],
    });
});

test('2D render models take field, forces and test-charge forces from one model', async ({ page }) => {
    const r = await page.evaluate(async () => {
        const s = await import('/solvers.js');
        const charges = [
            { x: 200, y: 150, charge: 2, id: 0 }, { x: 450, y: 180, charge: -1, id: 1 },
            { x: 300, y: 60, charge: 1, id: 2, isTestCharge: true },
        ];
        const out = {};
        for (const solver of ['electrostatic_line_charge', 'magnetostatic_wire']) {
            const m = s.makeRenderModel(solver, { field: null, charges, uniformX: 0.3, uniformY: -0.2, width: 668, height: 335 });
            const b = new Float64Array(2);
            m.fieldToBuffer(null, 320, 200, b);
            const f = m.model2D.fieldAt(320, 200);
            const force = m.forceOn(320, 200, 1.5);
            out[solver] = {
                hasModel: !!m.model2D, dimension: m.dimension,
                bufferMatches: Math.abs(b[0] - f.u) + Math.abs(b[1] - f.v),
                force, field: f, sources: m.model2D.sources.length,
                forces: m.forcesOnSources(),
            };
        }
        return out;
    });
    const e = r.electrostatic_line_charge, w = r.magnetostatic_wire;
    expect(e.hasModel && w.hasModel).toBe(true);
    expect(e.sources).toBe(2); // test charges don't generate field
    expect(e.bufferMatches).toBeLessThan(1e-12);
    expect(w.bufferMatches).toBeLessThan(1e-12);
    // electric force along E; magnetic force on a test current perpendicular to B (as in the 3D model)
    expect(e.force.u).toBeCloseTo(1.5 * e.field.u, 10);
    expect(e.force.v).toBeCloseTo(1.5 * e.field.v, 10);
    expect(w.force.u).toBeCloseTo(1.5 * w.field.v, 10);
    expect(w.force.v).toBeCloseTo(-1.5 * w.field.u, 10);
    // opposite charges attract; antiparallel currents repel
    expect(e.forces[0].u).toBeGreaterThan(0);
    expect(w.forces[0].u).toBeLessThan(0);
});

test('2D model normalisation matches the 3D field at a quarter of the canvas width (D3)', async ({ page }) => {
    const r = await page.evaluate(async () => {
        const s = await import('/solvers.js');
        const { FieldConfiguration } = await import('/maxwell/out/maxwell.js'); // initialised in beforeEach
        const W = 668, H = 335;
        const charges = [{ x: 100, y: 150, charge: 1, id: 0 }];
        const field = new FieldConfiguration(W, H, 64, 64);
        field.set_charges(charges);
        const m3 = s.makeRenderModel('electrostatic_direct', { field, charges, uniformX: 0, uniformY: 0, width: W, height: H });
        const m2 = s.makeRenderModel('electrostatic_line_charge', { field, charges, uniformX: 0, uniformY: 0, width: W, height: H });
        const b3 = new Float64Array(2), b2 = new Float64Array(2);
        m3.fieldToBuffer(field, 100 + W / 4, 150, b3);
        m2.fieldToBuffer(field, 100 + W / 4, 150, b2);
        return { e3: b3[0], e2: b2[0] };
    });
    expect(r.e2 / r.e3).toBeCloseTo(1, 3);
});

test('stream function and potential are consistent with the field', async ({ page }) => {
    const r = await page.evaluate(async () => {
        const { makeFieldModel } = await import('/fieldmodels.js');
        const charges = [{ x: 200, y: 150, charge: 2 }, { x: 450, y: 180, charge: -1 }];
        const e = makeFieldModel('electrostatic_line_charge', { charges, uniformX: 0.4, uniformY: 0.1, width: 668, height: 335 });
        const w = makeFieldModel('magnetostatic_wire', { charges, uniformX: 0.4, uniformY: 0.1, width: 668, height: 335 });
        const h = 1e-4, x = 330, y = 260;
        const d = (f, dx, dy) => (f(x + dx, y + dy) - f(x - dx, y - dy)) / (2 * h);
        const E = e.fieldAt(x, y), B = w.fieldAt(x, y);
        return {
            fromPsi: [2 * Math.PI * e.K * d(e.streamFunctionAt, 0, h), -2 * Math.PI * e.K * d(e.streamFunctionAt, h, 0)],
            fromPhi: [-d(e.potentialAt, h, 0), -d(e.potentialAt, 0, h)],
            E: [E.u, E.v],
            BdotGradA: B.u * d(w.potentialAt, h, 0) + B.v * d(w.potentialAt, 0, h),
            Bmag: Math.hypot(B.u, B.v),
        };
    });
    for (const k of [0, 1]) {
        expect(r.fromPhi[k]).toBeCloseTo(r.E[k], 4);
        expect(Math.abs(r.fromPsi[k] - r.E[k])).toBeLessThan(1e-3 * Math.hypot(...r.E)); // psi is unsoftened
    }
    expect(Math.abs(r.BdotGradA)).toBeLessThan(1e-6 * r.Bmag ** 2);
});

test('field lines per unit charge and the level offset: defaults, and set per embed', async ({ page }) => {
    const r = await page.evaluate(async () => {
        const wasm = await import('/maxwell/out/maxwell.js');
        const s = await import('/solvers.js');
        const { electricFieldLines } = await import('/fieldlines2d.js');
        const W = 668, H = 335;
        const charges = [{ x: 200, y: 167.5, charge: 2 }, { x: 468, y: 167.5, charge: -2 }];
        const ends = (lines) => lines.filter(l => l.from === 0).length;   // lines leaving the +2 charge
        const run = (solver, fieldLines) => {
            const field = wasm.FieldConfiguration.with_cell_size(W, H, 6);
            field.set_charges(charges);
            field.reset_fields();
            const m = s.makeRenderModel(solver, { field, charges, uniformX: 0, uniformY: 0, width: W, height: H, fieldLines }).model2D;
            const lines = electricFieldLines(m, { width: W, height: H });
            return { n: m.linesPerUnitCharge, f: m.levelOffset, ends: ends(lines), levels: [...new Set(lines.map(l => +l.level.toFixed(6)))].sort((a, b) => a - b) };
        };
        const out = {};
        for (const solver of ['electrostatic_line_charge', 'electrostatic_fourier'])
            out[solver] = [run(solver, {}), run(solver, { linesPerUnitCharge: 3, levelOffset: 0 }), run(solver, { linesPerUnitCharge: null })];
        return out;
    });
    for (const [solver, [byDefault, set, unset]] of Object.entries(r)) {
        expect(byDefault, solver).toMatchObject({ n: 5, f: 0.5, ends: 10 });
        // levels (k + 1/2) / 5
        for (const l of byDefault.levels) expect(Math.abs(5 * l - 0.5 - Math.round(5 * l - 0.5)), `${solver} level ${l}`).toBeLessThan(1e-5);
        expect(set, solver).toMatchObject({ n: 3, f: 0, ends: 6 });
        for (const l of set.levels) expect(Math.abs(3 * l - Math.round(3 * l)), `${solver} level ${l}`).toBeLessThan(1e-5);
        expect(unset, solver).toMatchObject({ n: 5, f: 0.5 });
    }
});

test('the field-line settings in an embed\'s state are used', async ({ page }) => {
    const { examples } = await import('./helpers.js');
    const base = { ...examples.find(e => e.id === 'lecture04-2').state, plotType: 'fieldline' };
    const render = async (state) => {
        await page.goto('/tests/harness.html?attrs=' + encodeURIComponent(JSON.stringify({ meme: JSON.stringify(state) })));
        const canvas = page.locator('meme-embed canvas').first();
        await expect(canvas).toBeVisible();
        await page.waitForTimeout(300);
        return canvas.evaluate(c => c.toDataURL());
    };
    const byDefault = await render(base);
    const explicit = await render({ ...base, linesPerUnitCharge: 5, fieldLineOffset: 0.5 });
    const original = await render({ ...base, linesPerUnitCharge: 4, fieldLineOffset: 0 });
    expect(explicit).toBe(byDefault);
    expect(original).not.toBe(byDefault);
});
