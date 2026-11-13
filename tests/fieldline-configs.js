// Source configurations used by the field-line rule tests. Pure data / seeded randomness; works
// in Node and in the browser.
//
// A configuration is {name, charges: [{x, y, charge, isTestCharge?}], uniformX, uniformY, width, height}
// in canvas CSS pixels, ready to pass to makeFieldModel(kind, config).

export const CANVAS_W = 668;
export const CANVAS_H = 335;

/** Small seeded PRNG (mulberry32); returns a function giving uniform numbers in [0, 1). */
export function rng(seed) {
    let a = seed >>> 0;
    return () => {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

/**
 * Configurations of the lecture examples (tests/fixtures/examples.json) with the given solver:
 * 'electrostatic_direct' (for the line-charge model) or 'magnetostatic_direct' (for the wire
 * model). Charges are stored in normalised coordinates (x relative to the canvas width,
 * y to the height); uniformElecX/Y are stored multiplied by the canvas height.
 */
export function lectureConfigs(examples, solver) {
    return examples.filter(e => e.state.solver === solver).map(e => ({
        name: e.id,
        charges: e.state.charges.map(c => ({
            x: c.x * CANVAS_W, y: c.y * CANVAS_H, charge: c.charge, isTestCharge: !!c.isTestCharge,
        })),
        uniformX: (e.state.uniformElecX || 0) / CANVAS_H,
        uniformY: (e.state.uniformElecY || 0) / CANVAS_H,
        width: CANVAS_W, height: CANVAS_H,
    }));
}

/**
 * `count` seeded random configurations of 3 or 4 charges of strength +-1..+-3, at least
 * `minSep` px apart and at least `margin` px from the edge of the canvas.
 */
export function randomConfigs(count = 20, { seed = 1234, minSep = 50, margin = 40, width = CANVAS_W, height = CANVAS_H } = {}) {
    const rand = rng(seed);
    const configs = [];
    for (let i = 0; i < count; i++) {
        const n = 3 + (i % 2); // alternate 3 and 4 charges
        const charges = [];
        let guard = 0;
        while (charges.length < n && guard++ < 10000) {
            const x = margin + rand() * (width - 2 * margin);
            const y = margin + rand() * (height - 2 * margin);
            if (charges.every(c => Math.hypot(c.x - x, c.y - y) >= minSep)) {
                const q = (1 + Math.floor(rand() * 3)) * (rand() < 0.5 ? -1 : 1);
                charges.push({ x: Math.round(x * 10) / 10, y: Math.round(y * 10) / 10, charge: q });
            }
        }
        configs.push({ name: `random-${i}`, charges, uniformX: 0, uniformY: 0, width, height });
    }
    return configs;
}

const cfg = (name, charges, width = CANVAS_W, height = CANVAS_H) =>
    ({ name, charges, uniformX: 0, uniformY: 0, width, height });

/**
 * Configurations for the drag test: `move` = {index, dx, dy} is the total displacement of one charge
 * (applied in 1 px steps). The moving charge stays at least 30 px from every other charge
 * throughout (checked by the test).
 */
export function dragConfigs() {
    return [
        { ...cfg('drag-triple', [{ x: 150, y: 167, charge: 2 }, { x: 500, y: 167, charge: -2 }, { x: 330, y: 70, charge: 1 }]),
          move: { index: 0, dx: 40, dy: 0 } },
        { ...cfg('drag-dipole-3', [{ x: 250, y: 100, charge: 3 }, { x: 450, y: 220, charge: -3 }]),
          move: { index: 0, dx: 0, dy: 40 } },
        { ...cfg('drag-quad', [{ x: 150, y: 80, charge: 1 }, { x: 500, y: 80, charge: -1 }, { x: 200, y: 260, charge: -2 }, { x: 450, y: 250, charge: 2 }]),
          move: { index: 3, dx: -28, dy: -28 } },
        { ...cfg('drag-heavy-neg', [{ x: 330, y: 167, charge: 3 }, { x: 500, y: 100, charge: -1 }, { x: 160, y: 250, charge: -1 }]),
          move: { index: 1, dx: -40, dy: 0 } },
        { ...cfg('drag-approach', [{ x: 170, y: 167, charge: 2 }, { x: 252, y: 167, charge: -2 }]),
          move: { index: 0, dx: 40, dy: 0 } }, // ends 42 px from the other charge
        { ...cfg('drag-like', [{ x: 200, y: 120, charge: 2 }, { x: 330, y: 200, charge: 2 }, { x: 450, y: 120, charge: -1 }]),
          move: { index: 1, dx: 40, dy: 0 } },
    ];
}

/** Neutral configurations on a 400 x 400 canvas for the rotation test (D5) */
export function rotationConfigs() {
    const S = 400;
    const c = (name, charges) => cfg(name, charges, S, S);
    return [
        c('rot-dipole-2', [{ x: 130, y: 200, charge: 2 }, { x: 270, y: 200, charge: -2 }]),
        c('rot-dipole-offaxis', [{ x: 100, y: 120, charge: 2 }, { x: 300, y: 240, charge: -2 }]),
        c('rot-dipole-1', [{ x: 150, y: 170, charge: 1 }, { x: 260, y: 250, charge: -1 }]),
        c('rot-two-plus-minus-2', [{ x: 120, y: 150, charge: 1 }, { x: 280, y: 150, charge: 1 }, { x: 200, y: 260, charge: -2 }]),
        c('rot-triple', [{ x: 110, y: 160, charge: 3 }, { x: 290, y: 120, charge: -1 }, { x: 230, y: 290, charge: -2 }]),
        c('rot-quadrupole', [{ x: 140, y: 140, charge: 1 }, { x: 260, y: 140, charge: -1 }, { x: 260, y: 260, charge: 1 }, { x: 140, y: 260, charge: -1 }]),
    ];
}

/** Rotate every charge of a configuration by `quarterTurns` x 90 degrees about the canvas centre. */
export function rotateConfig(config, quarterTurns = 1) {
    const cx = config.width / 2, cy = config.height / 2;
    const charges = config.charges.map(c => {
        let dx = c.x - cx, dy = c.y - cy;
        for (let i = 0; i < ((quarterTurns % 4) + 4) % 4; i++) [dx, dy] = [-dy, dx];
        return { ...c, x: cx + dx, y: cy + dy };
    });
    return { ...config, name: `${config.name}-rot${quarterTurns * 90}`, charges };
}

/** Horizontal axis pairs on the standard canvas */
export function axisConfigs() {
    return {
        unlike: cfg('axis-unlike', [{ x: 200, y: 167.5, charge: 2 }, { x: 468, y: 167.5, charge: -2 }]),
        like: cfg('axis-like', [{ x: 200, y: 167.5, charge: 2 }, { x: 468, y: 167.5, charge: 2 }]),
    };
}

export function totalCharge(config) {
    return config.charges.filter(c => !c.isTestCharge).reduce((t, c) => t + c.charge, 0);
}
