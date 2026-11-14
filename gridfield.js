// Field-line model for the grid solvers (electrostatic_fourier, dynamic).
// See https://github.com/apontzen/maxwell/issues/3 (stage 2b) and maxwell/src/streamfn.rs.
//
// The grid solvers are 2D, so their field lines are contours of a stream function psi. We use
//     psi = psi_sing + psi_rest + psi_uniform,
// psi_sing = sum q theta / 2 pi about the sources' true positions (multivalued; unwrapped by
// levelsets2d.js exactly as for the analytic line-charge model), psi_uniform the analytic uniform
// field term, and psi_rest a single-valued function computed by the WASM module on the corners of
// the grid's physical region (the canvas) and interpolated here (bicubic Catmull-Rom, so its
// gradient is continuous).
//
// The additive constant of psi_rest chooses which contours are drawn. For the dynamic solver the
// WASM module evolves it causally from frame to frame, so that lines ahead of the wavefront don't
// move when a charge is dragged (see maxwell/src/streamfn.rs).
//
// The model has the interface of fieldmodels.js's electric model as far as fieldlines2d.js needs
// it, plus `psiRest`. The world ends at the edge of the canvas (beyond it is the absorbing layer
// of the grid), so lines leaving the canvas end there ("infinity").

import { compute_grid_stream_function_rest, grid_physical_cells } from './maxwell/out/maxwell.js';
import { LINES_PER_UNIT_CHARGE, LEVEL_OFFSET } from './fieldmodels.js';

export const GRID_KIND = 'electrostatic_grid';

// Must match charge_normalization * cell area in maxwell/src/lib.rs (FieldConfiguration::new)
export const GRID_FLUX_PER_UNIT_CHARGE = 4000.0;

/** Bicubic (Catmull-Rom) interpolation of corner data (mx * my, x-major); mirrored (Neumann) beyond the edges. */
export class CornerInterpolant {
    constructor(data, mx, my, dx, dy) {
        this.data = data; this.mx = mx; this.my = my; this.dx = dx; this.dy = dy;
        this.val = 0; this.gx = 0; this.gy = 0;
        this.wx = new Float64Array(4); this.wy = new Float64Array(4);
        this.dwx = new Float64Array(4); this.dwy = new Float64Array(4);
    }
    static weights(t, w, dw) {
        const t2 = t * t, t3 = t2 * t;
        w[0] = 0.5 * (-t3 + 2 * t2 - t);   dw[0] = 0.5 * (-3 * t2 + 4 * t - 1);
        w[1] = 0.5 * (3 * t3 - 5 * t2 + 2); dw[1] = 0.5 * (9 * t2 - 10 * t);
        w[2] = 0.5 * (-3 * t3 + 4 * t2 + t); dw[2] = 0.5 * (-9 * t2 + 8 * t + 1);
        w[3] = 0.5 * (t3 - t2);             dw[3] = 0.5 * (3 * t2 - 2 * t);
    }
    /** Sets this.val, this.gx, this.gy at (x, y). */
    eval(x, y) {
        const mx = this.mx, my = this.my;
        let u = x / this.dx, v = y / this.dy;
        // clamp to the region (the tracer never needs values further out)
        u = Math.min(Math.max(u, 0), mx - 1); v = Math.min(Math.max(v, 0), my - 1);
        let i = Math.floor(u), j = Math.floor(v);
        if (i > mx - 2) i = mx - 2;
        if (j > my - 2) j = my - 2;
        CornerInterpolant.weights(u - i, this.wx, this.dwx);
        CornerInterpolant.weights(v - j, this.wy, this.dwy);
        const d = this.data;
        let f = 0, fx = 0, fy = 0;
        for (let a = 0; a < 4; a++) {
            let ia = i - 1 + a;
            if (ia < 0) ia = -1 - ia; else if (ia >= mx) ia = 2 * mx - 1 - ia;
            const row = ia * my;
            let s = 0, sy = 0;
            for (let b = 0; b < 4; b++) {
                let jb = j - 1 + b;
                if (jb < 0) jb = -1 - jb; else if (jb >= my) jb = 2 * my - 1 - jb;
                const z = d[row + jb];
                s += this.wy[b] * z; sy += this.dwy[b] * z;
            }
            f += this.wx[a] * s; fx += this.dwx[a] * s; fy += this.wx[a] * sy;
        }
        this.val = f; this.gx = fx / this.dx; this.gy = fy / this.dy;
    }
}

/**
 * @param {{field: FieldConfiguration, charges, uniformX, uniformY, width, height, linesPerUnitCharge,
 *          levelOffset}} config
 *   field must already hold the charges (and, for the dynamic solver, the current fields).
 *   linesPerUnitCharge, levelOffset: as for fieldmodels.js's electric model.
 */
export function makeGridFieldModel({ field, charges, uniformX = 0, uniformY = 0, width, height,
                                     linesPerUnitCharge = LINES_PER_UNIT_CHARGE, levelOffset = LEVEL_OFFSET }) {
    const sources = charges.filter(c => !c.isTestCharge);
    const [nx, ny] = grid_physical_cells(field);
    const xs = Float64Array.from(sources, c => c.x);
    const ys = Float64Array.from(sources, c => c.y);
    const qs = Float64Array.from(sources, c => c.charge);
    field.set_lines_per_unit_charge(linesPerUnitCharge);   // psi's constant is chosen modulo the line spacing
    const data = compute_grid_stream_function_rest(field, xs, ys, qs);
    const psiRest = new CornerInterpolant(data, nx + 1, ny + 1, width / nx, height / ny);
    const K = GRID_FLUX_PER_UNIT_CHARGE / (2 * Math.PI);
    // E from the stream function (the field the lines follow): E = 2 pi K rot(grad psi). Used for
    // the direction arrows and Gaussian-surface crossings, so that they agree with the lines.
    const fieldAt = (x, y) => {
        let u = uniformX, v = uniformY;
        for (const c of sources) {
            const dx = x - c.x, dy = y - c.y, f = K * c.charge / (dx * dx + dy * dy);
            u += f * dx; v += f * dy;
        }
        psiRest.eval(x, y);
        u += 2 * Math.PI * K * psiRest.gy;
        v -= 2 * Math.PI * K * psiRest.gx;
        return { u, v };
    };
    return {
        kind: GRID_KIND, isMagnetic: false, dimension: 2, sources,
        // E of a unit point charge at distance r is K / r (flux 2 pi K per unit charge)
        K,
        uniform: { x: uniformX, y: uniformY }, width, height,
        linesPerUnitCharge, levelOffset,
        fieldAt,
        fieldToBuffer(_f, x, y, buffer) {
            const { u, v } = fieldAt(x, y);
            buffer[0] = u; buffer[1] = v;
        },
        psiRest,
        // the world ends at the edge of the canvas
        region: { x0: 0, y0: 0, x1: width, y1: height },
    };
}
