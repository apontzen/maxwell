// Analytic 2D field models (line charges / straight wires), used for everything
// drawn in field-line mode and for the 2D solvers in any mode.
// See https://github.com/apontzen/maxwell/issues/3 (decisions D1-D6).
//
// Coordinates are CSS pixels of the canvas, as for the charges in ui.js.
//
// Every model has this interface:
//   kind                      'electrostatic_line_charge' | 'magnetostatic_wire'
//   isMagnetic                boolean
//   dimension                 2
//   sources                   the non-test charges/currents that generate the field
//   K                         strength constant: field of a unit source at distance r is K/r
//   uniform                   {x, y} uniform field (in field units, as stored by ui.js)
//   width, height             canvas size
//   fieldAt(x, y)             -> {u, v}  (E, or B for magnetic models)
//   fieldToBuffer(f, x, y, buffer)          same, in the signature of the WASM *_to_buffer functions
//                                           (first argument ignored), for draw.js's quiver code
//   fieldPerChargeToBuffer(f, x, y, buffer) per-source contributions (2 entries per source, no
//                                           uniform field), as compute_field_*_per_charge_direct_to_buffer
//   forceOn(x, y, q)          -> {u, v}  force on a test charge / test current q at (x, y)
//   forcesOnSources()         -> [{u, v}] force on each source from the other sources (excluding the
//                             uniform field, as the WASM compute_forces_* functions do)
//   potentialAt(x, y)         electric: the 2D potential phi; magnetic: A_z.
//                             For magnetic models the field lines are contours of A_z.
//   linesPerUnitCharge        electric only: field lines per unit charge (default LINES_PER_UNIT_CHARGE)
//   levelOffset               electric only: the drawn levels of psi are (k + levelOffset) / linesPerUnitCharge
//                             (default LEVEL_OFFSET). With 0, a +-q pair has a straight line running out to
//                             infinity along its axis on either side; with 1/2 it does not.
//   streamFunctionAt(x, y)    electric only: psi, normalised so that it increases by q on going once
//                             around a charge q in the direction of increasing atan2(y - y_i, x - x_i).
//                             Field lines are contours of psi at the levels (k + levelOffset) / linesPerUnitCharge.
//                             NB psi is multivalued: this returns the principal value, with each
//                             charge's branch cut along the ray theta = +-pi (pointing in -x from the
//                             charge). Callers tracing contours must unwrap.
//
// Conventions (mirroring the WASM 3D functions so that switching model doesn't flip anything):
//   electric:  E = K sum_i q_i (r - r_i) / |r - r_i|^2  + uniform
//   magnetic:  B = rotate(E_equivalent) with (u, v) -> (v, -u), where E_equivalent is the electric
//              field obtained by treating currents as charges (as compute_field_magnetostatic_direct)
//   force on a test current I in B = (u, v):  I * (v, -u)   (as compute_one_force_magnetostatic)
//   forces between currents: minus the equivalent electrostatic forces (as compute_forces_magnetostatic)
//
// Normalisation (decision D3): K = FIELD_SCALING_3D / r_ref with r_ref = width / 4, so that a
// 2D source and a 3D point source have the same field strength at distance r_ref.

export const FIELD_SCALING_3D = 20000.0; // must match FIELD_SCALING in maxwell/src/lib.rs
export const SOFTEN = 5.0;               // px^2; must match SOFTEN in maxwell/src/lib.rs
// Defaults for the electric field lines (see linesPerUnitCharge, levelOffset above; decision D4,
// revised: originally 4 and 0). Configurable per embed (ui.js: state.linesPerUnitCharge,
// state.fieldLineOffset).
export const LINES_PER_UNIT_CHARGE = 5;
export const LEVEL_OFFSET = 0.5;

export const LINE_CHARGE_KIND = 'electrostatic_line_charge';
export const WIRE_KIND = 'magnetostatic_wire';

function equivalentElectric(sources, K, ux, uy, x, y) {
    let u = ux, v = uy;
    for (const c of sources) {
        const dx = x - c.x, dy = y - c.y;
        const r2 = dx * dx + dy * dy + SOFTEN;
        u += K * c.charge * dx / r2;
        v += K * c.charge * dy / r2;
    }
    return { u, v };
}

function equivalentPotential(sources, K, ux, uy, width, height, x, y) {
    // phi = -K sum q ln(r / width) - uniform . (r - centre)
    let phi = -(ux * (x - width / 2) + uy * (y - height / 2));
    for (const c of sources) {
        const dx = x - c.x, dy = y - c.y;
        phi -= K * c.charge * 0.5 * Math.log((dx * dx + dy * dy + SOFTEN) / (width * width));
    }
    return phi;
}

/**
 * @param {string} kind  LINE_CHARGE_KIND or WIRE_KIND
 * @param {{charges: Array<{x:number,y:number,charge:number,isTestCharge?:boolean}>,
 *          uniformX?: number, uniformY?: number, width: number, height: number}} config
 */
export function makeFieldModel(kind, { charges, uniformX = 0, uniformY = 0, width, height,
                                       linesPerUnitCharge = LINES_PER_UNIT_CHARGE, levelOffset = LEVEL_OFFSET }) {
    if (kind !== LINE_CHARGE_KIND && kind !== WIRE_KIND)
        throw new Error(`Unknown 2D field model ${kind}`);
    const isMagnetic = kind === WIRE_KIND;
    const sources = charges.filter(c => !c.isTestCharge);
    const K = FIELD_SCALING_3D / (width / 4);
    const ux = uniformX, uy = uniformY;
    const xc = width / 2, yc = height / 2;

    const fieldAt = isMagnetic
        ? (x, y) => { const e = equivalentElectric(sources, K, ux, uy, x, y); return { u: e.v, v: -e.u }; }
        : (x, y) => equivalentElectric(sources, K, ux, uy, x, y);

    const model = {
        kind, isMagnetic, dimension: 2, sources, K, uniform: { x: ux, y: uy }, width, height,
        ...(!isMagnetic && { linesPerUnitCharge, levelOffset }),
        fieldAt,
        fieldToBuffer(_f, x, y, buffer) {
            const { u, v } = fieldAt(x, y);
            buffer[0] = u; buffer[1] = v;
        },
        fieldPerChargeToBuffer(_f, x, y, buffer) {
            sources.forEach((c, i) => {
                const e = equivalentElectric([c], K, 0, 0, x, y);
                if (isMagnetic) { buffer[2 * i] = e.v; buffer[2 * i + 1] = -e.u; }
                else { buffer[2 * i] = e.u; buffer[2 * i + 1] = e.v; }
            });
        },
        forceOn(x, y, q) {
            const { u, v } = fieldAt(x, y);
            return isMagnetic ? { u: v * q, v: -u * q } : { u: u * q, v: v * q };
        },
        forcesOnSources() {
            // Like compute_forces_electrostatic / _magnetostatic, this excludes the uniform field
            return sources.map(c => {
                const others = sources.filter(o => o !== c);
                const e = equivalentElectric(others, K, 0, 0, c.x, c.y);
                const s = isMagnetic ? -c.charge : c.charge;
                return { u: e.u * s, v: e.v * s };
            });
        },
        potentialAt(x, y) {
            const phi = equivalentPotential(sources, K, ux, uy, width, height, x, y);
            return isMagnetic ? -phi : phi;   // A_z = -phi_equivalent
        },
        streamFunctionAt: isMagnetic ? undefined : (x, y) => {
            let psi = (ux * (y - yc) - uy * (x - xc)) / (2 * Math.PI * K);
            for (const c of sources) psi += c.charge * Math.atan2(y - c.y, x - c.x) / (2 * Math.PI);
            return psi;
        },
    };
    return model;
}
