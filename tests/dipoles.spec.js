// Unit tests (Node) of the pairing of charges into rigid dipoles (dipoles.js), which happens when
// the "Dipoles" tickbox is ticked.

import { test, expect } from '@playwright/test';
import { pairDipoles, dropDanglingPartners, nextChargeId } from '../dipoles.js';

const q = (id, x, charge, extra = {}) => ({ id, x, y: 0, charge, ...extra });
/** The pairs, as sorted [id, id] */
const pairs = charges => charges.filter(c => c.dipoleWith !== undefined && c.id < c.dipoleWith).map(c => [c.id, c.dipoleWith]);

test.describe('pairDipoles', () => {
    test('each charge is paired with the nearest free charge of exactly -q', () => {
        const cs = [q(0, 0, 1), q(1, 50, -1), q(2, 20, -1), q(3, 100, 2), q(4, 130, -2)];
        expect(pairDipoles(cs)).toBe(2);
        expect(pairs(cs)).toEqual([[0, 2], [3, 4]]);
        // both ends refer to each other
        expect(cs[2].dipoleWith).toBe(0);
        expect(cs[1].dipoleWith).toBeUndefined();
    });

    test('greedily, in the order of the charges: a later charge may lose its nearest counterpart', () => {
        // +1 at 0 comes first and takes the only -1, although the +1 at 12 is nearer to it
        const cs = [q(0, 0, 1), q(1, 12, 1), q(2, 10, -1)];
        pairDipoles(cs);
        expect(pairs(cs)).toEqual([[0, 2]]);
        expect(cs[1].dipoleWith).toBeUndefined();
    });

    test('a charge without an exact counterpart stays free', () => {
        const cs = [q(0, 0, 2), q(1, 10, -1), q(2, 20, 1), q(3, 30, 0), q(4, 40, 0)];
        expect(pairDipoles(cs)).toBe(1);
        expect(pairs(cs)).toEqual([[1, 2]]);
        expect(cs[0].dipoleWith).toBeUndefined();
        // zero charges are not paired with each other
        expect(cs[3].dipoleWith).toBeUndefined();
    });

    test('charges already in dipoles keep them and are not considered', () => {
        const cs = [q(0, 0, 1, { dipoleWith: 1 }), q(1, 100, -1, { dipoleWith: 0 }), q(2, 5, -1), q(3, 95, 1)];
        expect(pairDipoles(cs)).toBe(1);
        expect(pairs(cs)).toEqual([[0, 1], [2, 3]]);
    });

    test('test charges pair only with test charges', () => {
        const cs = [q(0, 0, 1, { isTestCharge: true }), q(1, 10, -1), q(2, 50, -1, { isTestCharge: true }), q(3, 20, 1)];
        pairDipoles(cs);
        expect(pairs(cs)).toEqual([[0, 2], [1, 3]]);
    });

    test('distances are in 2D', () => {
        const cs = [{ id: 0, x: 0, y: 0, charge: 1 }, { id: 1, x: 30, y: 0, charge: -1 }, { id: 2, x: 0, y: 25, charge: -1 }];
        pairDipoles(cs);
        expect(pairs(cs)).toEqual([[0, 2]]);
    });

    test('broken references are dropped and duplicate or missing ids replaced before pairing', () => {
        const cs = [q(0, 0, 1, { dipoleWith: 7 }), q(0, 10, -1), { x: 20, y: 0, charge: 1 }, q(1, 30, -1)];
        pairDipoles(cs);
        const ids = cs.map(c => c.id);
        expect(new Set(ids).size).toBe(4);
        expect(pairs(cs).length).toBe(2);
        for (const c of cs) expect(cs.find(o => o.id === c.dipoleWith).dipoleWith).toBe(c.id);
    });
});

test('dropDanglingPartners frees the partner of a deleted charge; nextChargeId is unused', () => {
    const cs = [q(0, 0, 1, { dipoleWith: 1 }), q(1, 10, -1, { dipoleWith: 0 }), q(5, 20, 1)];
    const remaining = cs.filter(c => c.id !== 1);
    dropDanglingPartners(remaining);
    expect(remaining[0].dipoleWith).toBeUndefined();
    expect(nextChargeId(cs)).toBe(6);
    expect(nextChargeId([])).toBe(0);
});
