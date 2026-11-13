// Exact analytic field lines and equipotentials of a +q / -q pair of 2D line charges, used to
// validate the checkers (fieldline-checks.spec.js) and the rule runners without needing the real
// field-line implementation.
//
// For a +q charge at P and -q at M (field K q r / r^2 each):
//   * the field lines are circular arcs through both charges. The line leaving P at angle
//     alpha to the axis P->M is the circle through P and M whose tangent at P makes angle alpha
//     (measured with the same handedness as atan2) with the axis. With 4|q| lines at
//     alpha_k = 2 pi k / (4|q|) one line (k = 0) is the straight segment P->M and one (k = 2|q|, alpha =
//     pi) is a straight ray from P to infinity, matched by a ray from infinity into M.
//   * the equipotentials are the Apollonius circles |r - P| / |r - M| = lambda, exactly
//     orthogonal to every field line.

import { LINES_PER_UNIT_CHARGE } from '../fieldmodels.js';
import { TAU, normAngle } from './fieldline-checks.js';

const FAR = 1e5; // "infinity" for rays, in px

/**
 * @param {{x:number,y:number}} P  position of the positive charge
 * @param {{x:number,y:number}} M  position of the negative charge
 * @param {number} absq   |q|, so that there are 4|q| lines per charge
 * @param {{from?:number,to?:number,stepDeg?:number}} opts  indices of P and M in the source list
 * @returns lines {points, from, to, level} oriented from P to M (or to/from infinity)
 */
export function exactDipoleLines(P, M, absq, { from = 0, to = 1, stepDeg = 0.25 } = {}) {
    const n = Math.round(LINES_PER_UNIT_CHARGE * absq);
    const D = Math.hypot(M.x - P.x, M.y - P.y);
    const ux = (M.x - P.x) / D, uy = (M.y - P.y) / D;
    const nx = -uy, ny = ux;
    const lines = [];
    for (let k = 0; k < n; k++) {
        const alpha = (TAU * k) / n;
        const level = k / LINES_PER_UNIT_CHARGE;
        const tx = Math.cos(alpha) * ux + Math.sin(alpha) * nx;
        const ty = Math.cos(alpha) * uy + Math.sin(alpha) * ny;
        if (k === 0) {
            lines.push({ points: [[P.x, P.y], [M.x, M.y]], from, to, level });
        } else if (2 * k === n) {
            // along the axis, away from the other charge: to infinity from P, and from infinity into M
            lines.push({ points: [[P.x, P.y], [P.x + tx * FAR, P.y + ty * FAR]], from, to: 'infinity', level });
            lines.push({ points: [[M.x - tx * FAR, M.y - ty * FAR], [M.x, M.y]], from: 'infinity', to, level });
        } else {
            const h = -(D / 2) * Math.cos(alpha) / Math.sin(alpha);
            const cx = (P.x + M.x) / 2 + h * nx, cy = (P.y + M.y) / 2 + h * ny;
            const rho = Math.hypot(P.x - cx, P.y - cy);
            const phiP = Math.atan2(P.y - cy, P.x - cx), phiM = Math.atan2(M.y - cy, M.x - cx);
            // direction of travel at P is t; the circle is traversed counter-clockwise if the
            // cross product (P - C) x t is positive
            const s = (P.x - cx) * ty - (P.y - cy) * tx > 0 ? 1 : -1;
            const sweep = s * (s > 0 ? normAngle(phiM - phiP) : normAngle(phiP - phiM));
            const steps = Math.max(2, Math.ceil(Math.abs(sweep) / (stepDeg * Math.PI / 180)));
            const points = [];
            for (let i = 0; i <= steps; i++) {
                const phi = phiP + (sweep * i) / steps;
                points.push([cx + rho * Math.cos(phi), cy + rho * Math.sin(phi)]);
            }
            points[0] = [P.x, P.y];
            points[steps] = [M.x, M.y];
            lines.push({ points, from, to, level });
        }
    }
    return lines;
}

/**
 * Analytic attachment angle (chord angle at distance d, as defined by attachmentAngle) of the
 * exact dipole line with index k at the positive charge P.
 */
export function exactAttachmentAngleAtPlus(P, M, absq, k, d = 15) {
    const n = Math.round(LINES_PER_UNIT_CHARGE * absq);
    const D = Math.hypot(M.x - P.x, M.y - P.y);
    const axis = Math.atan2(M.y - P.y, M.x - P.x);
    const alpha = (TAU * k) / n;
    if (k === 0 || 2 * k === n) return normAngle(axis + alpha); // straight
    const rho = Math.abs(D / (2 * Math.sin(alpha)));
    const s = curvatureSign(P, M, alpha);
    const theta = 2 * Math.asin(Math.min(1, d / (2 * rho))); // arc angle subtended by the chord
    return normAngle(axis + alpha + s * theta / 2);
}

function curvatureSign(P, M, alpha) {
    // the circle centre is at (D/2, h) relative to P in axis coordinates, h = -(D/2) cot(alpha);
    // the travel direction t = (cos a, sin a); counter-clockwise if (P - C) x t > 0
    const D = Math.hypot(M.x - P.x, M.y - P.y);
    const h = -(D / 2) * Math.cos(alpha) / Math.sin(alpha);
    const cxr = D / 2, cyr = h;
    return (-cxr) * Math.sin(alpha) - (-cyr) * Math.cos(alpha) > 0 ? 1 : -1;
}

/**
 * Apollonius circle |r - P| / |r - M| = lambda as a closed polyline (or a long straight line for
 * lambda = 1). Returns {points, closed}.
 */
export function apolloniusCircle(P, M, lambda, { stepDeg = 1 } = {}) {
    const D = Math.hypot(M.x - P.x, M.y - P.y);
    const ux = (M.x - P.x) / D, uy = (M.y - P.y) / D;
    if (Math.abs(lambda - 1) < 1e-12) {
        const mx = (P.x + M.x) / 2, my = (P.y + M.y) / 2;
        return { points: [[mx - uy * FAR, my + ux * FAR], [mx + uy * FAR, my - ux * FAR]], closed: false };
    }
    const l2 = lambda * lambda;
    const sc = (-l2 * D) / (1 - l2); // centre, distance from P along the axis
    const rho = Math.sqrt(sc * sc + (l2 * D * D) / (1 - l2));
    const cx = P.x + sc * ux, cy = P.y + sc * uy;
    const n = Math.ceil(360 / stepDeg);
    const points = [];
    for (let i = 0; i < n; i++) {
        const phi = (TAU * i) / n;
        points.push([cx + rho * Math.cos(phi), cy + rho * Math.sin(phi)]);
    }
    points.push(points[0]);
    return { points, closed: true };
}
