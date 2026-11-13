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

import { LINES_PER_UNIT_CHARGE } from './fieldmodels.js';

export { LINES_PER_UNIT_CHARGE };

export function contourStep(model) {
    return model.K; // placeholder; stage 2 chooses the value
}

export function electricFieldLines(model, view) {
    return []; // stage 2
}

export function magneticFieldLines(model, view) {
    return []; // stage 2
}

export function potentialContours(model, view) {
    return []; // stage 2
}
