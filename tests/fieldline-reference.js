// A tiny EXACT reference implementation of the fieldlines2d.js contract, valid only for a neutral
// pair of line charges / wires (+q and -q, no uniform field). It exists so that the rule runners
// (fieldline-rules-lib.js) can be shown to pass on a correct implementation (see
// fieldline-checks.spec.js), i.e. that the rules are satisfiable and not over-strict.
//
// It is NOT the real implementation and is not used by the browser tests.

import { exactDipoleLines, apolloniusCircle } from './fieldline-exact.js';

const requirePair = model => {
    const s = model.sources;
    if (s.length !== 2 || Math.abs(s[0].charge + s[1].charge) > 1e-9 || model.uniform.x !== 0 || model.uniform.y !== 0)
        throw new Error('reference implementation only handles a neutral pair of sources');
    return s;
};

export function contourStep(model) {
    return model.K / 4;
}

export function electricFieldLines(model, view) {
    const s = requirePair(model);
    const iP = s[0].charge > 0 ? 0 : 1, iM = 1 - iP;
    return exactDipoleLines(s[iP], s[iM], Math.abs(s[iP].charge),
        { from: iP, to: iM, lpu: model.linesPerUnitCharge, levelOffset: model.levelOffset });
}

function circles(model, sign) {
    const s = requirePair(model);
    const step = contourStep(model);
    // potential phi = -K q ln(rP/rM) for source 0 = P; A_z = +K q ln(rP/rM)
    const q = s[0].charge, out = [];
    const kmax = Math.floor((3.5 * model.K * Math.abs(q)) / step);
    for (let k = -kmax; k <= kmax; k++) {
        const level = k * step;
        const lambda = Math.exp((-sign * level) / (model.K * q));
        const c = apolloniusCircle(s[0], s[1], lambda);
        out.push({ points: c.points, closed: c.closed, level });
    }
    return out;
}

export const potentialContours = (model, view) => circles(model, +1);
export const magneticFieldLines = (model, view) => circles(model, -1);
