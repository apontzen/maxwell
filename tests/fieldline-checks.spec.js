// Unit tests for the field-line checker functions (fieldline-checks.js) and the rule runners
// (fieldline-rules-lib.js), using synthetic but EXACT data. These run in Node, need no browser
// and no field-line implementation, and must pass at all times: they establish that the checkers
// themselves can be trusted before they are used on the real implementation
// (fieldline-rules.spec.js).
//
// The exact data are the field lines of a +q / -q pair of 2D line charges: circular arcs through
// both charges (see fieldline-exact.js), and their orthogonal Apollonius-circle equipotentials.

import { test, expect } from '@playwright/test';
import * as C from './fieldline-checks.js';
import { exactDipoleLines, exactAttachmentAngleAtPlus, apolloniusCircle } from './fieldline-exact.js';
import * as configs from './fieldline-configs.js';
import * as rules from './fieldline-rules-lib.js';
import * as reference from './fieldline-reference.js';
import * as models from '../fieldmodels.js';

const { DEG } = C;
const near = (a, b, tol, msg) => expect(Math.abs(a - b), msg).toBeLessThanOrEqual(tol);

/** Charges +q at P and -q at M as a source list */
const pair = (P, M, q) => [{ x: P.x, y: P.y, charge: q }, { x: M.x, y: M.y, charge: -q }];

const PAIRS = [
    { name: 'horizontal q=1', P: { x: 200, y: 167 }, M: { x: 468, y: 167 }, q: 1 },
    { name: 'horizontal q=2', P: { x: 200, y: 167 }, M: { x: 468, y: 167 }, q: 2 },
    { name: 'tilted q=3', P: { x: 250, y: 90 }, M: { x: 400, y: 230 }, q: 3 },
    { name: 'close q=2 (60 px)', P: { x: 300, y: 170 }, M: { x: 330, y: 222 }, q: 2 },
    { name: 'negative direction q=4', P: { x: 500, y: 250 }, M: { x: 180, y: 100 }, q: 4 },
];

test.describe('exact data generators', () => {
    test('exact dipole lines have the right number of lines and labels', () => {
        for (const { name, P, M, q } of PAIRS) {
            const lines = exactDipoleLines(P, M, q);
            // 4q lines leave P: 4q-2 arcs, the axis segment and a ray to infinity; plus a ray from infinity into M
            expect(lines.length, name).toBe(4 * q + 1);
            expect(lines.filter(l => l.from === 0 && l.to === 1).length, name).toBe(4 * q - 1);
            expect(lines.filter(l => l.to === 'infinity').length, name).toBe(1);
            expect(lines.filter(l => l.from === 'infinity').length, name).toBe(1);
        }
    });

    test('exact dipole arcs are tangent to the exact field (sanity check of the generator itself)', () => {
        for (const { name, P, M, q } of PAIRS) {
            const exactField = (x, y) => {
                const a = { u: x - P.x, v: y - P.y }, b = { u: x - M.x, v: y - M.y };
                const ra = a.u * a.u + a.v * a.v, rb = b.u * b.u + b.v * b.v;
                return { u: q * (a.u / ra - b.u / rb), v: q * (a.v / ra - b.v / rb) };
            };
            for (const l of exactDipoleLines(P, M, q)) {
                const finite = l.from === 0 && l.to === 1;
                if (!finite) continue;
                const r = C.maxTangentMisalignment(l.points, exactField, { oriented: true, sources: pair(P, M, q), minDist: 0.5 });
                expect(r.max / DEG, `${name} level ${l.level}`).toBeLessThan(0.2);
            }
        }
    });
});

test.describe('attachment counts', () => {
    for (const { name, P, M, q } of PAIRS) {
        test(`4|q| line ends per source: ${name}`, () => {
            const lines = exactDipoleLines(P, M, q);
            expect(C.countAttachedEnds(lines, pair(P, M, q), 1)).toEqual([4 * q, 4 * q]);
        });
    }

    test('an end within 1 px counts, an end 1.5 px away does not', () => {
        const S = [{ x: 100, y: 100, charge: 1 }];
        const mk = (dx) => [{ points: [[100 + dx, 100], [150, 100]], from: 0, to: 'infinity', level: 0 }];
        expect(C.countAttachedEnds(mk(0.9), S, 1)).toEqual([1]);
        expect(C.countAttachedEnds(mk(1.5), S, 1)).toEqual([0]);
    });

    test('both ends of a line are counted, and a missing line is detected', () => {
        const { P, M, q } = PAIRS[1];
        const lines = exactDipoleLines(P, M, q);
        const dropped = lines.filter((_, i) => i !== 3);
        expect(C.countAttachedEnds(dropped, pair(P, M, q), 1)).toEqual([4 * q - 1, 4 * q - 1]);
        // a duplicate line adds an end to each charge
        expect(C.countAttachedEnds([...lines, lines[2]], pair(P, M, q), 1)).toEqual([4 * q + 1, 4 * q + 1]);
    });
});

test.describe('attachment angles', () => {
    for (const { name, P, M, q } of PAIRS) {
        test(`angles at the positive charge match the analytic chord angles: ${name}`, () => {
            const lines = exactDipoleLines(P, M, q);
            const got = C.attachmentAngles(lines, pair(P, M, q), { d: 15, tol: 1 });
            const n = 4 * q;
            const expected = [];
            for (let k = 0; k < n; k++) expected.push(exactAttachmentAngleAtPlus(P, M, q, k, 15));
            expected.sort((a, b) => a - b);
            expect(got[0].length).toBe(n);
            const r = C.matchAngleSets(got[0], expected);
            expect(r.maxDiff / DEG, `max difference ${r.maxDiff / DEG} deg`).toBeLessThan(0.01);
        });

        test(`as d -> 0 the attachment angle tends to alpha_k = 2 pi k / 4|q| (the tangent): ${name}`, () => {
            const lines = exactDipoleLines(P, M, q, { stepDeg: 0.02 });
            const axis = Math.atan2(M.y - P.y, M.x - P.x);
            const got = C.attachmentAngles(lines, pair(P, M, q), { d: 0.05, tol: 1 })[0];
            const expected = Array.from({ length: 4 * q }, (_, k) => C.normAngle(axis + (2 * Math.PI * k) / (4 * q)));
            const r = C.matchAngleSets(got, expected);
            expect(r.maxDiff / DEG, `max difference ${r.maxDiff / DEG} deg`).toBeLessThan(0.1);
        });
    }

    test('angles at the negative charge are the mirror image of those at the positive charge', () => {
        // horizontal pair: reflection in the perpendicular bisector sends angle t to pi - t
        const { P, M, q } = PAIRS[1];
        const A = C.attachmentAngles(exactDipoleLines(P, M, q), pair(P, M, q));
        const mirrored = A[0].map(t => C.normAngle(Math.PI - t));
        const r = C.matchAngleSets(mirrored, A[1]);
        expect(r.maxDiff / DEG).toBeLessThan(0.01);
    });

    test('a line along the axis has attachment angle exactly towards the other charge', () => {
        const { P, M, q } = PAIRS[2];
        const angles = C.attachmentAngles(exactDipoleLines(P, M, q), pair(P, M, q))[0];
        const axis = Math.atan2(M.y - P.y, M.x - P.x);
        expect(Math.min(...angles.map(a => Math.abs(C.angleDiff(a, axis)))) / DEG).toBeLessThan(1e-9);
        // and one pointing exactly away from it (the ray to infinity)
        expect(Math.min(...angles.map(a => Math.abs(C.angleDiff(a, axis + Math.PI)))) / DEG).toBeLessThan(1e-9);
    });

    test('attachment angle of a line that ends at the source is measured walking away from that end', () => {
        // straight line arriving at the source from the left: outward direction is pi
        const S = { x: 100, y: 100, charge: -1 };
        const pts = [[0, 100], [50, 100], [100, 100]];
        expect(C.attachmentAngle(pts, 'end', S, 15) / DEG).toBeCloseTo(180, 9);
        // and one leaving upwards (y decreasing): angle -90 = 270
        expect(C.attachmentAngle([[100, 100], [100, 0]], 'start', S, 15) / DEG).toBeCloseTo(270, 9);
    });

    test('a line shorter than d falls back to the far end', () => {
        const S = { x: 0, y: 0, charge: 1 };
        expect(C.attachmentAngle([[0, 0], [3, 3]], 'start', S, 15) / DEG).toBeCloseTo(45, 9);
    });
});

test.describe('angle set matching', () => {
    const base = [10, 100, 130, 250, 355].map(d => d * DEG);

    test('identical sets match exactly, in any cyclic order', () => {
        const rotated = [...base.slice(2), ...base.slice(0, 2)];
        expect(C.matchAngleSets(base, rotated).maxDiff).toBeLessThan(1e-12);
    });

    test('a rotation by an offset is recovered', () => {
        const B = base.map(a => C.normAngle(a + 90 * DEG));
        expect(C.matchAngleSets(base, B, 90 * DEG).maxDiff).toBeLessThan(1e-12);
        expect(C.matchAngleSets(base, B, 0).maxDiff / DEG).toBeGreaterThan(10);
    });

    test('wrap-around at 0 / 360 degrees does not matter', () => {
        const A = [359, 1, 180].map(d => d * DEG), B = [1.5, 358.5, 180].map(d => d * DEG);
        expect(C.matchAngleSets(A, B).maxDiff / DEG).toBeCloseTo(0.5, 9);
    });

    test('a single perturbed angle is reported with the right size', () => {
        const B = base.map((a, i) => i === 3 ? a + 2 * DEG : a);
        expect(C.matchAngleSets(base, B).maxDiff / DEG).toBeCloseTo(2, 9);
    });

    test('sets of different sizes do not match', () => {
        const r = C.matchAngleSets(base, base.slice(1));
        expect(r.ok).toBe(false);
        expect(r.maxDiff).toBe(Infinity);
    });
});

test.describe('signed crossings with a polygon (Gauss)', () => {
    const square = (cx, cy, h) => [[cx - h, cy - h], [cx + h, cy - h], [cx + h, cy + h], [cx - h, cy + h]];
    const squareCW = (cx, cy, h) => square(cx, cy, h).reverse();

    test('simple geometry: leaving is +1, entering is -1, independent of polygon orientation', () => {
        for (const mk of [square, squareCW]) {
            const poly = mk(0, 0, 10);
            expect(C.signedCrossings([[0, 0], [30, 0]], poly)).toEqual({ out: 1, in: 0, net: 1 });
            expect(C.signedCrossings([[30, 0], [0, 0]], poly)).toEqual({ out: 0, in: 1, net: -1 });
            expect(C.signedCrossings([[-30, 0], [30, 0]], poly)).toEqual({ out: 1, in: 1, net: 0 });   // one segment crossing twice
            expect(C.signedCrossings([[-30, 20], [30, 20]], poly).net).toBe(0);                       // misses
            expect(C.signedCrossings([[0, 0], [5, 5]], poly).net).toBe(0);                            // stays inside
        }
    });

    test('a crossing exactly through a polygon vertex is counted once', () => {
        const poly = square(0, 0, 10);
        expect(C.signedCrossings([[0, 0], [20, 20]], poly).net).toBe(1);        // through the corner (10, 10)
        expect(C.signedCrossings([[0, 0], [10, 10], [20, 20]], poly).net).toBe(1);   // vertex of the polyline on the corner
        expect(C.signedCrossings([[0, 0], [10, 0], [20, 0]], poly).net).toBe(1);     // polyline vertex on an edge
    });

    test('flux of exact dipole lines = 4 x enclosed charge for polygons around +, -, both, neither', () => {
        for (const { name, P, M, q } of PAIRS) {
            const lines = exactDipoleLines(P, M, q);
            const mid = { x: (P.x + M.x) / 2, y: (P.y + M.y) / 2 };
            const r = Math.hypot(P.x - M.x, P.y - M.y);
            const around = (c, h) => square(c.x, c.y, h);
            // squares of half-size 9 (charges at the centre, so >= 9 px from the boundary and nothing else nearby)
            const cases = [
                [around(P, Math.min(9, r / 4)), 4 * q, 'around +'],
                [around(M, Math.min(9, r / 4)), -4 * q, 'around -'],
                [around(mid, r * 2), 0, 'around both'],
                [around({ x: mid.x + 3 * r, y: mid.y + 3 * r }, 20), 0, 'around nothing'],
                [around({ x: P.x - 2 * r, y: P.y - 2 * r }, 30), 0, 'far away'],
            ];
            for (const [poly, expected, what] of cases)
                expect(C.netCrossings(lines, poly).net, `${name}, ${what}`).toBe(expected);
            // half-plane-like huge polygon that contains + only (exactly one line family leaves)
            const tall = [[P.x - 5, P.y - 1e4], [P.x + 5, P.y - 1e4], [P.x + 5, P.y + 1e4], [P.x - 5, P.y + 1e4]];
            const inside = C.enclosedCharge(pair(P, M, q), tall);
            expect(C.netCrossings(lines, tall).net, `${name}, tall thin box`).toBe(4 * inside);
        }
    });

    test('rays to infinity are crossed once, lines between charges inside net to zero', () => {
        const { P, M, q } = PAIRS[0];
        const lines = exactDipoleLines(P, M, q);
        const poly = square((P.x + M.x) / 2, P.y, 300);
        // all finite lines start and end inside: crossings cancel; the two rays: + goes out (+1), - ray comes in (-1)
        const c = C.netCrossings(lines, poly);
        expect(c.net).toBe(0);
        expect(c.out).toBe(c.in);
    });

    test('random polygons against random exact dipoles: Gauss counts are exact', () => {
        const rand = configs.rng(7);
        for (let trial = 0; trial < 40; trial++) {
            const q = 1 + Math.floor(rand() * 4);
            const P = { x: 100 + rand() * 450, y: 60 + rand() * 220 };
            const M = { x: 100 + rand() * 450, y: 60 + rand() * 220 };
            if (Math.hypot(P.x - M.x, P.y - M.y) < 50) continue;
            const S = pair(P, M, q);
            const lines = exactDipoleLines(P, M, q);
            const cx = rand() * 668, cy = rand() * 335;
            const poly = [];
            const n = 5 + Math.floor(rand() * 5);
            for (let i = 0; i < n; i++) {
                const t = ((i + 0.5 * (rand() - 0.5)) / n) * 2 * Math.PI, r = 30 + rand() * 250;
                poly.push([cx + r * Math.cos(t), cy + r * Math.sin(t)]);
            }
            if (S.some(s => C.pointPolygonDistance([s.x, s.y], poly) < 15)) continue;
            expect(C.netCrossings(lines, poly).net, JSON.stringify({ trial, P, M, q, poly })).toBe(4 * C.enclosedCharge(S, poly));
        }
    });
});

test.describe('angles between polylines', () => {
    test('straight lines cross at the right angle (folded to 0..90 degrees, direction ignored)', () => {
        for (const a of [10, 45, 90, 135, 170]) {
            const t = a * DEG;
            const A = [[-10, 0], [10, 0]], B = [[-10 * Math.cos(t), -10 * Math.sin(t)], [10 * Math.cos(t), 10 * Math.sin(t)]];
            const folded = Math.min(a, 180 - a);
            const x = C.crossingAngles(A, B);
            expect(x.length).toBe(1);
            expect(x[0] / DEG).toBeCloseTo(folded, 9);
            expect(C.crossingAngles(A, [...B].reverse())[0] / DEG).toBeCloseTo(folded, 9);
        }
    });

    test('non-crossing polylines have no crossings; a zigzag crossing three times is found three times', () => {
        expect(C.crossingAngles([[0, 0], [10, 0]], [[0, 5], [10, 5]])).toEqual([]);
        expect(C.crossingAngles([[0, 0], [10, 0]], [[1, -1], [2, 1], [3, -1], [4, 1]]).length).toBe(3);
    });

    test('exact field lines and Apollonius equipotentials of a dipole are perpendicular where they cross', () => {
        for (const { name, P, M, q } of PAIRS.slice(0, 4)) {
            const lines = exactDipoleLines(P, M, q, { stepDeg: 0.1 });
            const contours = [];
            for (const lam of [0.05, 0.2, 0.5, 0.8, 1.25, 2, 5, 20]) contours.push(apolloniusCircle(P, M, lam, { stepDeg: 0.1 }));
            contours.push(apolloniusCircle(P, M, 1));
            const rect = { x0: -1500, y0: -1500, x1: 1500, y1: 1500 };
            const crossings = C.setCrossings(contours, lines, rect);
            expect(crossings.length, name).toBeGreaterThan(20);
            const worst = Math.max(...crossings.map(c => 90 - c.angle / DEG));
            expect(worst, `${name}: worst deviation from 90 deg`).toBeLessThan(0.2);
        }
    });

    test('setCrossings agrees with brute force polylineCrossings', () => {
        const { P, M, q } = PAIRS[2];
        const lines = exactDipoleLines(P, M, q, { stepDeg: 2 }).filter(l => l.points.length > 2);
        const contours = [0.3, 0.7, 1.4, 3].map(l => apolloniusCircle(P, M, l, { stepDeg: 3 }));
        let brute = 0;
        for (const a of contours) for (const b of lines) brute += C.polylineCrossings(a.points, b.points).length;
        expect(C.setCrossings(contours, lines, { x0: -1500, y0: -1500, x1: 1500, y1: 1500 }).length).toBe(brute);
    });

    test('a field line crosses a deliberately skewed contour at the wrong angle', () => {
        const A = [[0, -10], [0, 10]], B = [[-10, -10 * Math.tan(10 * DEG)], [10, 10 * Math.tan(10 * DEG)]]; // 10 deg off horizontal
        expect(90 - C.crossingAngles(A, B)[0] / DEG).toBeCloseTo(10, 9);
    });
});

test.describe('tangent misalignment', () => {
    // uniform field in direction (1, 1)
    const uniform = () => ({ u: 1, v: 1 });

    test('a line along the field has zero misalignment, one against it has 180 degrees (oriented) or 0 (unoriented)', () => {
        const along = [[0, 0], [5, 5], [10, 10]];
        expect(C.maxTangentMisalignment(along, uniform).max).toBeLessThan(1e-12);
        const against = [...along].reverse();
        expect(C.maxTangentMisalignment(against, uniform).max / DEG).toBeCloseTo(180, 9);
        expect(C.maxTangentMisalignment(against, uniform, { oriented: false }).max).toBeLessThan(1e-12);
    });

    test('a kink is reported with its size and location; sources and weak-field regions are skipped', () => {
        const line = [[0, 0], [5, 5], [10, 5]]; // second segment is 45 deg off
        const r = C.maxTangentMisalignment(line, uniform);
        expect(r.max / DEG).toBeCloseTo(45, 9);
        expect(r.at).toEqual([7.5, 5]);
        expect(C.maxTangentMisalignment(line, uniform, { sources: [{ x: 7.5, y: 5 }], minDist: 1 }).max / DEG).toBeCloseTo(0, 9);
        expect(C.maxTangentMisalignment(line, uniform, { minField: 10 }).n).toBe(0);
        expect(C.maxTangentMisalignment(line, uniform, { include: (x) => x < 6 }).max).toBeLessThan(1e-12);
    });

    test('the unsoftened field of a +-q pair is tangent to the exact arcs; the model field (softened) is within 3 degrees', () => {
        // fieldmodels.js uses 1/(r^2 + 5) in place of 1/r^2; away from the sources (> 5 px) this
        // should change the direction by well under the 3 degrees allowed by the rule tests
        for (const { name, P, M, q } of PAIRS) {
            const m = models.makeFieldModel(models.LINE_CHARGE_KIND, { charges: pair(P, M, q), width: 668, height: 335 });
            let worst = 0;
            for (const l of exactDipoleLines(P, M, q)) {
                const r = C.maxTangentMisalignment(l.points, m.fieldAt, { sources: m.sources, minDist: 5, minField: m.K / 400, include: (x, y) => Math.hypot(x - 334, y - 167) < 2 * Math.hypot(668, 335) });
                worst = Math.max(worst, r.max / DEG);
            }
            expect(worst, `${name}: worst misalignment with model.fieldAt`).toBeLessThan(3);
        }
    });

    test('a perturbed exact line is detected', () => {
        const { P, M, q } = PAIRS[1];
        const m = models.makeFieldModel(models.LINE_CHARGE_KIND, { charges: pair(P, M, q), width: 668, height: 335 });
        const l = exactDipoleLines(P, M, q).find(l => l.points.length > 100);
        const bad = l.points.map(([x, y], i) => [x + (i > 200 && i < 215 ? 3 : 0), y]);
        expect(C.maxTangentMisalignment(bad, m.fieldAt, { sources: m.sources }).max / DEG).toBeGreaterThan(3);
    });
});

test.describe('endpoint check', () => {
    const view = { width: 668, height: 335 };
    const { P, M, q } = PAIRS[1];
    const S = pair(P, M, q);

    test('exact lines pass', () => {
        expect(C.checkEndpoints(exactDipoleLines(P, M, q), S, view)).toEqual([]);
    });

    test('wrong labels, short rays and ends off the source are all reported', () => {
        const lines = exactDipoleLines(P, M, q);
        const swapped = lines.map((l, i) => i === 3 ? { ...l, from: 1 } : l);
        expect(C.checkEndpoints(swapped, S, view).length).toBeGreaterThan(0);
        const shortRay = lines.map(l => l.to === 'infinity' ? { ...l, points: [l.points[0], [l.points[0][0] - 500, l.points[0][1]]] } : l);
        expect(C.checkEndpoints(shortRay, S, view).join('\n')).toMatch(/infinity.*from the view centre/);
        const off = lines.map((l, i) => i === 2 ? { ...l, points: l.points.map((p, j) => j === 0 ? [p[0] + 3, p[1]] : p) } : l);
        expect(C.checkEndpoints(off, S, view).join('\n')).toMatch(/px from that source/);
        const mislabelled = lines.map((l, i) => i === 2 ? { ...l, to: 'infinity' } : l);
        expect(C.checkEndpoints(mislabelled, S, view).join('\n')).toMatch(/labelled 'infinity' but lies on source/);
    });
});

test.describe('comparing lines between views', () => {
    const { P, M, q } = PAIRS[2];
    const rect = { x0: 0, y0: 0, x1: 668, y1: 335 };

    test('different sampling of the same lines agree; a shifted line does not', () => {
        const a = exactDipoleLines(P, M, q, { stepDeg: 0.3 });
        const b = exactDipoleLines(P, M, q, { stepDeg: 0.37 });
        expect(C.compareLinesInRect(a, b, rect, { tol: 0.5 })).toEqual([]);
        const shifted = b.map((l, i) => i === 3 ? { ...l, points: l.points.map(([x, y], j) => j > 0 && j < l.points.length - 1 ? [x + 1, y] : [x, y]) } : l);
        expect(C.compareLinesInRect(a, shifted, rect, { tol: 0.5 }).length).toBeGreaterThan(0);
    });

    test('a line missing from one set is reported, and lines wholly outside the rect are ignored', () => {
        const a = exactDipoleLines(P, M, q), b = a.filter((_, i) => i !== 4);
        expect(C.compareLinesInRect(a, b, rect).join('\n')).toMatch(/no second line/);
        const far = { points: [[5000, 5000], [5100, 5100]], from: 'infinity', to: 'infinity', level: 7 };
        expect(C.compareLinesInRect([...a, far], a, rect)).toEqual([]);
    });
});

// ---------------------------------------------------------------------------
// The rule runners must hold for an exact implementation (so that they are satisfiable and the
// tolerances are not over-strict) and must fail for deliberately broken ones.
// ---------------------------------------------------------------------------

test.describe('rule runners against the exact reference (neutral pairs only)', () => {
    const env = { models, lines2d: reference };
    const pairCfg = (name, P, M, q, w = 668, h = 335) => ({ name, charges: [{ ...P, charge: q }, { ...M, charge: -q }], uniformX: 0, uniformY: 0, width: w, height: h });
    const PAIR_CFGS = PAIRS.map(p => pairCfg(p.name, p.P, p.M, p.q));
    const all = configs.dragConfigs();
    const dipoleDrags = all.filter(c => c.charges.length === 2);
    const dipoleRots = configs.rotationConfigs().filter(c => c.charges.length === 2);

    test('exact counts', () => expect(rules.ruleExactCounts(env, PAIR_CFGS)).toEqual([]));
    test('endpoints', () => expect(rules.ruleEndpoints(env, PAIR_CFGS)).toEqual([]));
    test('tangency (against the softened model.fieldAt)', () => expect(rules.ruleTangency(env, PAIR_CFGS)).toEqual([]));
    test('viewport independence', () => {
        expect(rules.ruleViewportIndependence(env, PAIR_CFGS)).toEqual([]);
        expect(rules.ruleViewportIndependenceContours(env, PAIR_CFGS, PAIR_CFGS)).toEqual([]);
    });
    test('drag stability', () => expect(rules.ruleDragStability(env, dipoleDrags)).toEqual([]));
    test('rotation invariance', () => expect(rules.ruleRotationInvariance(env, dipoleRots)).toEqual([]));
    test('axis line (unlike pair)', () => {
        const ax = configs.axisConfigs();
        const problems = rules.ruleAxisLine({ ...env, lines2d: reference }, { unlike: ax.unlike, like: { ...ax.unlike, name: 'skip' } });
        // the "like" half of the rule cannot be satisfied by a dipole-only reference: only the unlike half is asserted
        expect(problems.filter(p => p.startsWith('axis-unlike'))).toEqual([]);
    });
    test('contours perpendicular to field lines', () => expect(rules.rulePerpendicularContours(env, PAIR_CFGS)).toEqual([]));
    test("Gauss's law", () => expect(rules.ruleGauss(env, PAIR_CFGS)).toEqual([]));
    test('magnetic contours', () => expect(rules.ruleMagnetic(env, PAIR_CFGS)).toEqual([]));
    test('contour step scaling', () => expect(rules.ruleContourStepScaling(env, PAIR_CFGS)).toEqual([]));
    test('time budget', () => expect(rules.ruleTimeBudget(env, PAIR_CFGS.slice(0, 2))).toEqual([]));

    test('broken implementations are caught', () => {
        const base = { ...reference };
        const dropLine = { ...base, electricFieldLines: (m, v) => base.electricFieldLines(m, v).filter((_, i) => i !== 3) };
        expect(rules.ruleExactCounts({ models, lines2d: dropLine }, PAIR_CFGS).length).toBeGreaterThan(0);
        expect(rules.ruleGauss({ models, lines2d: dropLine }, PAIR_CFGS).length).toBeGreaterThan(0);
        const clipped = { ...base, electricFieldLines: (m, v) => base.electricFieldLines(m, v).map(l => ({ ...l, points: l.points.filter(([x, y]) => Math.abs(x - v.width / 2) < v.width && Math.abs(y - v.height / 2) < v.height) })) };
        expect(rules.ruleEndpoints({ models, lines2d: clipped }, PAIR_CFGS).length).toBeGreaterThan(0);
        const viewDependent = { ...base, electricFieldLines: (m, v) => base.electricFieldLines(m, v).map(l => ({ ...l, points: l.points.map(([x, y]) => [x + (v.width > 700 ? 1 : 0), y]) })) };
        expect(rules.ruleViewportIndependence({ models, lines2d: viewDependent }, PAIR_CFGS).length).toBeGreaterThan(0);
        const reversed = { ...base, electricFieldLines: (m, v) => base.electricFieldLines(m, v).map(l => ({ ...l, points: [...l.points].reverse(), from: l.to, to: l.from })) };
        expect(rules.ruleTangency({ models, lines2d: reversed }, PAIR_CFGS).length).toBeGreaterThan(0);
        expect(rules.ruleEndpoints({ models, lines2d: reversed }, PAIR_CFGS).length).toBeGreaterThan(0);
        const open = { ...base, magneticFieldLines: (m, v) => base.magneticFieldLines(m, v).map(l => ({ ...l, closed: false })) };
        // a neutral pair may have open contours only if they extend far outside the view; circles marked open do not
        expect(rules.ruleMagnetic({ models, lines2d: open }, [pairCfg('open-circles', PAIRS[0].P, PAIRS[0].M, 1)]).length).toBeGreaterThan(0);
        const slow = { ...base, electricFieldLines: (m, v) => { const t = Date.now(); while (Date.now() - t < 40); return base.electricFieldLines(m, v); } };
        expect(rules.ruleTimeBudget({ models, lines2d: slow }, PAIR_CFGS.slice(0, 1), { calls: 3 }).length).toBeGreaterThan(0);
        const stub = { electricFieldLines: () => [], potentialContours: () => [], magneticFieldLines: () => [], contourStep: m => m.K };
        for (const [name, f] of Object.entries({
            counts: () => rules.ruleExactCounts({ models, lines2d: stub }, PAIR_CFGS),
            endpoints: () => rules.ruleEndpoints({ models, lines2d: stub }, PAIR_CFGS),
            tangency: () => rules.ruleTangency({ models, lines2d: stub }, PAIR_CFGS),
            viewport: () => rules.ruleViewportIndependence({ models, lines2d: stub }, PAIR_CFGS),
            drag: () => rules.ruleDragStability({ models, lines2d: stub }, dipoleDrags),
            rotation: () => rules.ruleRotationInvariance({ models, lines2d: stub }, dipoleRots),
            perpendicular: () => rules.rulePerpendicularContours({ models, lines2d: stub }, PAIR_CFGS),
            gauss: () => rules.ruleGauss({ models, lines2d: stub }, PAIR_CFGS),
            magnetic: () => rules.ruleMagnetic({ models, lines2d: stub }, PAIR_CFGS),
            time: () => rules.ruleTimeBudget({ models, lines2d: stub }, PAIR_CFGS.slice(0, 1), { calls: 2 }),
            axis: () => rules.ruleAxisLine({ models, lines2d: stub }, configs.axisConfigs()),
            contours: () => rules.ruleViewportIndependenceContours({ models, lines2d: stub }, PAIR_CFGS, PAIR_CFGS),
        }))
            expect(f().length, `${name} rule must fail on the stub`).toBeGreaterThan(0);
    });
});

test.describe('configurations', () => {
    test('random configurations are reproducible, 3-4 charges of +-1..3, at least 50 px apart', () => {
        const a = configs.randomConfigs(20), b = configs.randomConfigs(20);
        expect(a).toEqual(b);
        expect(a.length).toBe(20);
        for (const c of a) {
            expect([3, 4]).toContain(c.charges.length);
            for (const q of c.charges) { expect(Math.abs(q.charge)).toBeGreaterThanOrEqual(1); expect(Math.abs(q.charge)).toBeLessThanOrEqual(3); expect(Number.isInteger(q.charge)).toBe(true); }
            for (let i = 0; i < c.charges.length; i++) for (let j = i + 1; j < c.charges.length; j++)
                expect(Math.hypot(c.charges[i].x - c.charges[j].x, c.charges[i].y - c.charges[j].y)).toBeGreaterThanOrEqual(50);
        }
    });

    test('drag configurations keep the moving charge at least 30 px from the others', () => {
        for (const c of configs.dragConfigs()) {
            const { index, dx, dy } = c.move;
            for (let s = 0; s <= 40; s++) {
                const x = c.charges[index].x + (dx * s) / 40, y = c.charges[index].y + (dy * s) / 40;
                for (let j = 0; j < c.charges.length; j++) if (j !== index)
                    expect(Math.hypot(c.charges[j].x - x, c.charges[j].y - y), `${c.name} step ${s}`).toBeGreaterThanOrEqual(30);
            }
        }
    });

    test('rotation configurations are neutral and stay on the 400 px canvas when rotated', () => {
        for (const c of configs.rotationConfigs()) {
            expect(configs.totalCharge(c)).toBe(0);
            for (const k of [1, 2, 3]) for (const q of configs.rotateConfig(c, k).charges) {
                expect(q.x).toBeGreaterThan(20); expect(q.x).toBeLessThan(380); expect(q.y).toBeGreaterThan(20); expect(q.y).toBeLessThan(380);
            }
        }
        // rotating four times is the identity
        const c = configs.rotationConfigs()[1];
        const r4 = configs.rotateConfig(c, 4);
        c.charges.forEach((q, i) => { near(q.x, r4.charges[i].x, 1e-9); near(q.y, r4.charges[i].y, 1e-9); });
    });

    test('lecture configurations convert normalised coordinates and the uniform field as documented', async () => {
        const fs = await import('node:fs');
        const ex = JSON.parse(fs.readFileSync(new URL('./fixtures/examples.json', import.meta.url), 'utf8'));
        const el = configs.lectureConfigs(ex, 'electrostatic_direct'), mag = configs.lectureConfigs(ex, 'magnetostatic_direct');
        expect(el.length).toBe(ex.filter(e => e.state.solver === 'electrostatic_direct').length);
        expect(mag.length).toBe(ex.filter(e => e.state.solver === 'magnetostatic_direct').length);
        const l32 = el.find(c => c.name === 'lecture03-2');
        near(l32.uniformX, 400 / 335, 1e-12);
        expect(l32.charges.every(c => c.isTestCharge)).toBe(true);
        const l41 = el.find(c => c.name === 'lecture04-1');
        near(l41.charges[0].x, 0.3 * 668, 1e-9); near(l41.charges[0].y, 0.5 * 335, 1e-9);
        // test charges are not sources of the model
        const m = models.makeFieldModel(models.WIRE_KIND, mag.find(c => c.name === 'lecture09-2'));
        expect(m.sources.length).toBe(1);
    });
});
