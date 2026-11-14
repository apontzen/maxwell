// Registry of solvers (physical models), and adapters giving every solver the same
// interface for drawing. See https://github.com/apontzen/maxwell/issues/3 (D1, D6).
//
// The solver stored in the state/JSON is the user's choice. In field-line mode a 3D
// solver is replaced by its 2D counterpart (effectiveSolver), and *everything* drawn
// then comes from that 2D model: lines, potential contours, forces, flux arrows.

import {
    compute_field_electrostatic_direct_to_buffer, compute_field_magnetostatic_direct_to_buffer,
    compute_field_electrostatic_per_charge_direct_to_buffer, compute_field_magnetostatic_per_charge_direct_to_buffer,
    compute_electric_field_dynamic_to_buffer, compute_forces_electrostatic, compute_forces_magnetostatic,
    compute_one_force_electrostatic, compute_one_force_magnetostatic,
} from './maxwell/out/maxwell.js';
import { makeFieldModel, LINE_CHARGE_KIND, WIRE_KIND } from './fieldmodels.js';
import { makeGridFieldModel } from './gridfield.js';

export const SOLVERS = {
    electrostatic_direct: {
        label: 'Point charges [3D, in-plane]', type: '3d', magnetic: false,
        fieldlineCounterpart: 'electrostatic_line_charge',
        potential: true, forces: true, perCharge: true, fieldlines: true,
    },
    electrostatic_line_charge: {
        label: 'Line charges [2D]', type: '2d', magnetic: false,
        potential: true, forces: true, perCharge: true, fieldlines: true,
    },
    electrostatic_fourier: {
        label: 'Line charges [2D, grid solver]', type: 'grid', magnetic: false,
        potential: false, forces: false, perCharge: false, fieldlines: true,
    },
    dynamic: {
        label: 'Dynamic (Maxwell) [2D]', type: 'grid', magnetic: false, dynamic: true,
        potential: false, forces: false, perCharge: false, fieldlines: true,
    },
    magnetostatic_moving_charge: {
        label: 'Moving charges [3D, in-plane]', type: '3d', magnetic: true,
        fieldlineCounterpart: 'magnetostatic_wire',
        potential: false, forces: true, perCharge: true, fieldlines: true,
    },
    magnetostatic_wire: {
        label: 'Straight wires [2D]', type: '2d', magnetic: true,
        potential: false, forces: true, perCharge: true, fieldlines: true,
    },
};

/** Deprecated names still accepted in JSON (D1). */
export const SOLVER_ALIASES = {
    magnetostatic_direct: 'magnetostatic_moving_charge',
};

export function canonicalSolver(name) {
    const canonical = SOLVER_ALIASES[name] ?? name;
    if (!(canonical in SOLVERS)) {
        console.error(`Unknown solver type ${name}`);
        return 'electrostatic_direct';
    }
    return canonical;
}

export function isMagnetic(name) {
    return SOLVERS[canonicalSolver(name)].magnetic;
}

/** The solver actually used for drawing, given the user's choice and the plot type. */
export function effectiveSolver(name, plotType) {
    const s = canonicalSolver(name);
    if (plotType === 'fieldline' && SOLVERS[s].fieldlineCounterpart)
        return SOLVERS[s].fieldlineCounterpart;
    return s;
}

/**
 * A uniform drawing interface over all solvers.
 *
 * @param {string} solver  an effective (canonical) solver name
 * @param {object} opts    {field: FieldConfiguration (WASM), charges (incl. test charges),
 *                          uniformX, uniformY, width, height, perCharge}
 * @returns {{solver, info, isMagnetic, dimension, perCharge, field, model2D,
 *            fieldToBuffer(f, x, y, buf), fieldPerChargeToBuffer(f, x, y, buf),
 *            forceOn(x, y, q), forcesOnSources()}}
 *   fieldToBuffer / fieldPerChargeToBuffer take the WASM-style (field, x, y, buffer) arguments;
 *   model2D is the analytic model from fieldmodels.js for 2D solvers, otherwise null.
 */
export function makeRenderModel(solver, { field, charges, uniformX, uniformY, width, height, perCharge }) {
    const info = SOLVERS[solver];
    const base = {
        solver, info, isMagnetic: info.magnetic, field, model2D: null,
        perCharge: !!(perCharge && info.perCharge),
    };

    if (info.type === '2d') {
        const model2D = makeFieldModel(info.magnetic ? WIRE_KIND : LINE_CHARGE_KIND,
            { charges, uniformX, uniformY, width, height });
        return {
            ...base, dimension: 2, model2D,
            fieldToBuffer: model2D.fieldToBuffer,
            fieldPerChargeToBuffer: model2D.fieldPerChargeToBuffer,
            forceOn: model2D.forceOn,
            forcesOnSources: model2D.forcesOnSources,
        };
    }

    if (info.type === '3d') {
        const m = info.magnetic;
        return {
            ...base, dimension: 3,
            fieldToBuffer: m ? compute_field_magnetostatic_direct_to_buffer : compute_field_electrostatic_direct_to_buffer,
            fieldPerChargeToBuffer: m ? compute_field_magnetostatic_per_charge_direct_to_buffer
                : compute_field_electrostatic_per_charge_direct_to_buffer,
            forceOn: (x, y, q) => (m ? compute_one_force_magnetostatic : compute_one_force_electrostatic)(field, x, y, q),
            forcesOnSources: () => (m ? compute_forces_magnetostatic : compute_forces_electrostatic)(field),
        };
    }

    // grid solvers (Fourier / dynamic)
    return {
        ...base, dimension: 2, perCharge: false,
        // built lazily: only field-line mode needs it
        get model2D() {
            if (!this._grid) {
                this._grid = makeGridFieldModel({ field, charges, uniformX, uniformY, width, height });
                // the dynamic solver's lines evolve from frame to frame: keep their arrows in place
                this._grid.stickyArrows = !!info.dynamic;
            }
            return this._grid;
        },
        fieldToBuffer: compute_electric_field_dynamic_to_buffer,
        fieldPerChargeToBuffer: compute_electric_field_dynamic_to_buffer,
        forceOn: () => ({ u: 0, v: 0 }),
        forcesOnSources: () => null,
    };
}
