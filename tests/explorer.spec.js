// Tests for the stand-alone explorer (index.html), in which charges can be added,
// edited, deleted and the solver changed.

import {
    test, expect, waitForRender, canvasPixels, countDifferentPixels, colourAt,
    classifyColour, advance, installPausedClock,
} from './helpers.js';

const CANVAS = '#vectorFieldCanvas';

/** A starting state, in the same JSON format used by the embeds and "Copy"/"Paste" */
const TWO_CHARGES = {
    charges: [
        { x: 0.3, y: 0.5, charge: 2, id: 0 },
        { x: 0.7, y: 0.5, charge: -1, id: 1 },
    ],
    solver: 'electrostatic_direct', showPotential: false, plotType: 'quiver', perCharge: false,
};

const DIPOLES = {
    charges: [
        { x: 0.4, y: 0.3, charge: 1, id: 0, isTestCharge: true, dipoleWith: 1 },
        { x: 0.6, y: 0.3, charge: -1, id: 1, isTestCharge: true, dipoleWith: 0 },
        { x: 0.4, y: 0.7, charge: 2, id: 2, dipoleWith: 3 },
        { x: 0.6, y: 0.7, charge: -2, id: 3, dipoleWith: 2 },
    ],
    solver: 'electrostatic_direct', showPotential: false, plotType: 'fieldline', perCharge: false,
    dipoleMode: true,
};

const SOLVERS = ['electrostatic_direct', 'electrostatic_line_charge', 'electrostatic_fourier', 'dynamic',
    'magnetostatic_moving_charge', 'magnetostatic_wire'];

/** Which optional controls the explorer shows for each solver */
const CONTROLS_FOR_SOLVER = {
    electrostatic_direct: { potential: true, fieldlines: true, forces: true, perCharge: true },
    electrostatic_line_charge: { potential: true, fieldlines: true, forces: true, perCharge: true },
    electrostatic_fourier: { potential: false, fieldlines: true, forces: false, perCharge: false },
    dynamic: { potential: false, fieldlines: true, forces: false, perCharge: false },
    magnetostatic_moving_charge: { potential: false, fieldlines: true, forces: true, perCharge: true },
    magnetostatic_wire: { potential: false, fieldlines: true, forces: true, perCharge: true },
};

/** The 2D model each 3D model switches to when field lines are shown (issue #3, D6) */
const FIELDLINE_COUNTERPART = {
    electrostatic_direct: 'electrostatic_line_charge',
    magnetostatic_moving_charge: 'magnetostatic_wire',
};

const CONTROL_IDS = {
    potential: '#show-potential-control', fieldlines: '#fieldlines-control',
    forces: '#forces-control', perCharge: '#per-charge-control',
};

async function snapshot(page) {
    const dataUrl = await page.locator(CANVAS).evaluate(c => c.toDataURL('image/png'));
    return Buffer.from(dataUrl.split(',')[1], 'base64');
}

const pixels = page => canvasPixels(page.locator(CANVAS));

/** Load the explorer, optionally with a saved state in localStorage */
async function openExplorer(page, state = null) {
    if (state) {
        await page.addInitScript(s => {
            if (!localStorage.getItem('maxwell_state')) localStorage.setItem('maxwell_state', s);
        }, JSON.stringify(state));
    }
    await page.goto('/index.html');
    await waitForRender(page, 1, { selector: CANVAS, allowBlank: !state });
}

async function canvasPoint(page, { x, y }) {
    const m = await page.locator(CANVAS).evaluate(c => {
        const r = c.getBoundingClientRect();
        return { left: r.left + c.clientLeft, top: r.top + c.clientTop, width: r.width, height: r.height };
    });
    return { x: Math.round(m.left + x * m.width), y: Math.round(m.top + y * m.height) };
}

async function drag(page, from, to, steps = 10) {
    const a = await canvasPoint(page, from), b = await canvasPoint(page, to);
    await page.mouse.move(a.x, a.y);
    await page.mouse.down();
    for (let i = 1; i <= steps; i++)
        await page.mouse.move(Math.round(a.x + (b.x - a.x) * i / steps), Math.round(a.y + (b.y - a.y) * i / steps));
    await page.mouse.up();
    await advance(page);
}

const colourAtPoint = (page, p) => colourAt(page.locator(CANVAS), p);

test.describe('explorer (index.html)', () => {

    test('starts empty and "Add Charge" adds a positive charge in the centre', async ({ page }) => {
        await openExplorer(page);
        await expect(page.locator('#solver')).toHaveValue('electrostatic_direct');
        await page.locator('#addPositiveCharge').click();
        await advance(page);
        expect(classifyColour(await colourAtPoint(page, { x: 0.5, y: 0.5 }))).toBe('red');
        expect(await snapshot(page)).toMatchSnapshot(['add-charge.png']);
    });

    test('double-click adds a charge whose properties can be edited, then deletes it', async ({ page }) => {
        await openExplorer(page, TWO_CHARGES);
        const initial = await pixels(page);
        const p = await canvasPoint(page, { x: 0.5, y: 0.25 });
        await page.mouse.dblclick(p.x, p.y);
        await advance(page);

        const props = page.locator('.charge-properties');
        await expect(props).toBeVisible();
        await expect(props.locator('.chargeValue')).toHaveText('1');
        expect(classifyColour(await colourAtPoint(page, { x: 0.5, y: 0.25 }))).toBe('red');
        expect(await snapshot(page)).toMatchSnapshot(['dblclick-add.png']);

        // Change its strength with the slider
        await props.locator('input.charge').fill('-3');
        await advance(page);
        await expect(props.locator('.chargeValue')).toHaveText('-3');
        expect(classifyColour(await colourAtPoint(page, { x: 0.5, y: 0.25 }))).toBe('blue');
        expect(await snapshot(page)).toMatchSnapshot(['edit-strength.png']);

        // Make it a test charge (drawn in a pale colour; generates no field)
        await props.locator('input.testCharge').check();
        await advance(page);
        expect(classifyColour(await colourAtPoint(page, { x: 0.5, y: 0.25 }))).toBe('other');
        expect(await snapshot(page)).toMatchSnapshot(['test-charge.png']);

        // Delete it
        await props.locator('.deleteCharge').click();
        await advance(page);
        await expect(props).toBeHidden();
        expect(countDifferentPixels(initial, await pixels(page))).toBe(0);
    });

    test('double-clicking an existing charge deletes it', async ({ page }) => {
        await openExplorer(page, TWO_CHARGES);
        const p = await canvasPoint(page, TWO_CHARGES.charges[1]);
        await page.mouse.dblclick(p.x, p.y);
        await advance(page);
        expect(classifyColour(await colourAtPoint(page, TWO_CHARGES.charges[1]))).not.toBe('blue');
        expect(await snapshot(page)).toMatchSnapshot(['dblclick-delete.png']);
    });

    test('a short click on a charge selects it without moving it', async ({ page }) => {
        await openExplorer(page, TWO_CHARGES);
        const p = await canvasPoint(page, TWO_CHARGES.charges[0]);
        await page.mouse.click(p.x, p.y);
        await advance(page);
        await expect(page.locator('.charge-properties')).toBeVisible();
        await expect(page.locator('.charge-properties .chargeValue')).toHaveText('2');
        // clicking elsewhere deselects
        const q = await canvasPoint(page, { x: 0.05, y: 0.05 });
        await page.mouse.click(q.x, q.y);
        await advance(page);
        await expect(page.locator('.charge-properties')).toBeHidden();
    });

    test('dragging a charge moves it', async ({ page }) => {
        await openExplorer(page, TWO_CHARGES);
        await drag(page, TWO_CHARGES.charges[0], { x: 0.4, y: 0.2 });
        expect(classifyColour(await colourAtPoint(page, { x: 0.4, y: 0.2 }))).toBe('red');
        expect(classifyColour(await colourAtPoint(page, TWO_CHARGES.charges[0]))).not.toBe('red');
        await expect(page.locator('.charge-properties')).toBeHidden();
        expect(await snapshot(page)).toMatchSnapshot(['drag.png']);
    });

    test('"Clear" removes all charges', async ({ page }) => {
        await openExplorer(page, TWO_CHARGES);
        await page.locator('#clearCharges').click();
        await advance(page);
        for (const c of TWO_CHARGES.charges)
            expect(classifyColour(await colourAtPoint(page, c))).not.toMatch(/red|blue/);
        expect(await snapshot(page)).toMatchSnapshot(['clear.png']);
    });

    test('state is restored from localStorage after a reload', async ({ page }) => {
        await openExplorer(page, TWO_CHARGES);
        await drag(page, TWO_CHARGES.charges[1], { x: 0.6, y: 0.8 });
        await page.locator('#potential').check();
        await advance(page);
        const before = await pixels(page);
        await page.reload();
        await waitForRender(page, 1, { selector: CANVAS });
        await expect(page.locator('#potential')).toBeChecked();
        expect(countDifferentPixels(before, await pixels(page))).toBe(0);
    });

    for (const solver of SOLVERS) {
        test(`solver "${solver}" shows the right controls and renders`, async ({ page }) => {
            await openExplorer(page, TWO_CHARGES);
            await page.locator('#solver').selectOption(solver);
            await advance(page, solver === 'dynamic' ? 1000 : 100);
            for (const [key, visible] of Object.entries(CONTROLS_FOR_SOLVER[solver])) {
                if (visible) await expect(page.locator(CONTROL_IDS[key]), key).toBeVisible();
                else await expect(page.locator(CONTROL_IDS[key]), key).toBeHidden();
            }
            expect(await snapshot(page)).toMatchSnapshot([`solver-${solver}.png`]);
        });
    }

    test('button labels say "Current" for the magnetostatic solvers', async ({ page }) => {
        await openExplorer(page, TWO_CHARGES);
        await page.locator('#solver').selectOption('magnetostatic_wire');
        await advance(page);
        await expect(page.locator('#addPositiveCharge')).toHaveText('Add Current', { timeout: 1000 });
        await page.locator('#solver').selectOption('magnetostatic_moving_charge');
        await advance(page);
        await expect(page.locator('#addPositiveCharge')).toHaveText('Add Current', { timeout: 1000 });
        await page.locator('#solver').selectOption('electrostatic_direct');
        await advance(page);
        await expect(page.locator('#addPositiveCharge')).toHaveText('Add Charge', { timeout: 1000 });
    });

    const TOGGLES = [
        { id: '#potential', solver: 'electrostatic_direct', state: TWO_CHARGES },
        { id: '#fieldlines', solver: 'electrostatic_direct', state: TWO_CHARGES },
        { id: '#forces', solver: 'electrostatic_direct', state: TWO_CHARGES },
        { id: '#per-charge', solver: 'electrostatic_direct', state: TWO_CHARGES },
        { id: '#potential', solver: 'electrostatic_line_charge', state: TWO_CHARGES },
        { id: '#fieldlines', solver: 'electrostatic_line_charge', state: TWO_CHARGES },
        { id: '#fieldlines', solver: 'magnetostatic_moving_charge', state: TWO_CHARGES },
        { id: '#forces', solver: 'magnetostatic_moving_charge', state: TWO_CHARGES },
        { id: '#per-charge', solver: 'magnetostatic_moving_charge', state: TWO_CHARGES },
        { id: '#fieldlines', solver: 'magnetostatic_wire', state: TWO_CHARGES },
        { id: '#forces', solver: 'magnetostatic_wire', state: TWO_CHARGES },
        { id: '#per-charge', solver: 'magnetostatic_wire', state: TWO_CHARGES },
        { id: '#dipole', solver: 'electrostatic_direct', state: DIPOLES },
    ];
    for (const { id, solver, state } of TOGGLES) {
        test(`toggle ${id} with solver "${solver}"`, async ({ page }) => {
            await openExplorer(page, { ...state, solver });
            const box = page.locator(id);
            const wasChecked = await box.isChecked();
            const before = await pixels(page);
            await box.click();
            await advance(page);
            expect(await box.isChecked()).toBe(!wasChecked);
            expect(await snapshot(page)).toMatchSnapshot([`toggle-${id.slice(1)}-${solver}.png`]);
            expect(countDifferentPixels(before, await pixels(page))).toBeGreaterThan(0);
            await box.click();
            await advance(page);
            expect(countDifferentPixels(before, await pixels(page))).toBe(0);
        });
    }

    for (const [solver3D, solver2D] of Object.entries(FIELDLINE_COUNTERPART)) {
        test(`ticking field lines visibly switches "${solver3D}" to "${solver2D}" and back`, async ({ page, context }) => {
            await context.grantPermissions(['clipboard-read', 'clipboard-write']);
            await openExplorer(page, { ...TWO_CHARGES, solver: solver3D });
            const dropdown = page.locator('#solver');
            await expect(dropdown).toHaveValue(solver3D);

            await page.locator('#fieldlines').check();
            await advance(page);
            await expect(dropdown).toHaveValue(solver2D);
            const viaSwitch = await pixels(page);

            // The saved/copied state records the user's choice plus the tickbox, not the derived model
            await page.locator('#copyJson').click();
            await advance(page);
            const copied = JSON.parse(await page.evaluate(() => navigator.clipboard.readText()));
            expect(copied.solver).toBe(solver3D);
            expect(copied.plotType).toBe('fieldline');

            await page.locator('#fieldlines').uncheck();
            await advance(page);
            await expect(dropdown).toHaveValue(solver3D);

            // Field lines drawn via the switch are identical to choosing the 2D model directly
            await dropdown.selectOption(solver2D);
            await page.locator('#fieldlines').check();
            await advance(page);
            expect(countDifferentPixels(viaSwitch, await pixels(page))).toBe(0);
            // ...and having explicitly chosen the 2D model, unticking keeps it
            await page.locator('#fieldlines').uncheck();
            await advance(page);
            await expect(dropdown).toHaveValue(solver2D);
        });
    }

    test('the deprecated "magnetostatic_direct" solver name means moving charges', async ({ page, browser }) => {
        await openExplorer(page, { ...TWO_CHARGES, solver: 'magnetostatic_direct' });
        await expect(page.locator('#solver')).toHaveValue('magnetostatic_moving_charge');
        const viaAlias = await pixels(page);

        const other = await browser.newContext({ viewport: page.viewportSize(), deviceScaleFactor: 1 });
        const page2 = await other.newPage();
        await installPausedClock(page2);
        await openExplorer(page2, { ...TWO_CHARGES, solver: 'magnetostatic_moving_charge' });
        expect(countDifferentPixels(viaAlias, await pixels(page2))).toBe(0);
        await other.close();
    });

    test('rigid dipoles rotate and translate as a unit', async ({ page }) => {
        await openExplorer(page, DIPOLES);
        // drag the end of the lower dipole upwards: rotates it about its centre
        await drag(page, DIPOLES.charges[2], { x: 0.4, y: 0.55 });
        expect(await snapshot(page)).toMatchSnapshot(['dipole-rotate.png']);
        // drag the midpoint of the upper (test) dipole: translates it
        await drag(page, { x: 0.5, y: 0.3 }, { x: 0.5, y: 0.15 });
        expect(await snapshot(page)).toMatchSnapshot(['dipole-translate.png']);
    });

    test('"Copy" puts the state on the clipboard and "Paste" restores it', async ({ page, context }) => {
        await context.grantPermissions(['clipboard-read', 'clipboard-write']);
        await openExplorer(page, TWO_CHARGES);
        await drag(page, TWO_CHARGES.charges[0], { x: 0.2, y: 0.3 });
        const dragged = await pixels(page);

        await page.locator('#copyJson').click();
        await advance(page);
        const copied = JSON.parse(await page.evaluate(() => navigator.clipboard.readText()));
        expect(copied.solver).toBe('electrostatic_direct');
        expect(copied.charges).toHaveLength(2);
        expect(copied.charges[0].x).toBeCloseTo(0.2, 1);
        expect(copied.charges[0].y).toBeCloseTo(0.3, 1);

        await page.locator('#clearCharges').click();
        await advance(page);
        await page.locator('#pasteJson').click();
        await advance(page);
        await expect.poll(async () => countDifferentPixels(dragged, await pixels(page))).toBe(0);
    });
});
