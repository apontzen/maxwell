// Pairing charges into rigid dipoles (ui.js: ticking the "Dipoles" / "Rigid dipole" tickbox).
//
// A dipole is two charges that refer to each other by id: a.dipoleWith === b.id and
// b.dipoleWith === a.id (this is also how a state's JSON specifies dipoles).

/**
 * Pair up the charges (in canvas coordinates; modified in place) greedily into dipoles: going
 * through the charges in order, each charge q not yet in a dipole is paired with the nearest charge
 * of exactly -q that is not in a dipole either (and is a test charge if and only if it is). A charge
 * with no such counterpart stays free. Charges already in dipoles are left as they are. Returns the
 * number of new dipoles.
 */
export function pairDipoles(charges) {
    ensureUniqueIds(charges);
    dropDanglingPartners(charges);
    const byId = new Map(charges.map(c => [c.id, c]));
    const inDipole = c => c.dipoleWith !== undefined && byId.get(c.dipoleWith)?.dipoleWith === c.id;
    let made = 0;
    for (const c of charges) {
        if (inDipole(c) || c.charge === 0) continue;
        let best = null, bestD2 = Infinity;
        for (const o of charges) {
            if (o === c || inDipole(o) || o.charge !== -c.charge || !!o.isTestCharge !== !!c.isTestCharge) continue;
            const d2 = (o.x - c.x) ** 2 + (o.y - c.y) ** 2;
            if (d2 < bestD2) { bestD2 = d2; best = o; }
        }
        if (best) {
            c.dipoleWith = best.id;
            best.dipoleWith = c.id;
            made++;
        }
    }
    return made;
}

/** Remove dipoleWith from charges whose partner is missing or doesn't refer back (e.g. after a deletion). */
export function dropDanglingPartners(charges) {
    const byId = new Map(charges.map(c => [c.id, c]));
    for (const c of charges)
        if (c.dipoleWith !== undefined && byId.get(c.dipoleWith)?.dipoleWith !== c.id) delete c.dipoleWith;
}

/** An id for a new charge: one more than the largest in use. */
export function nextChargeId(charges) {
    return charges.reduce((m, c) => (Number.isFinite(c.id) ? Math.max(m, c.id) : m), -1) + 1;
}

/** Give charges without an id, or sharing one with an earlier charge, a fresh id (none of these can be in a dipole). */
function ensureUniqueIds(charges) {
    const seen = new Set();
    let next = nextChargeId(charges);
    for (const c of charges) {
        if (c.id === undefined || seen.has(c.id)) {
            delete c.dipoleWith;
            c.id = next++;
        }
        seen.add(c.id);
    }
}
