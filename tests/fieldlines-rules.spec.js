// Rule tests for the 2D field-line geometry (fieldlines2d.js), run on the real implementation in
// the browser. The rules themselves live in fieldline-rules-lib.js and are built from the checkers
// in fieldline-checks.js, which fieldline-checks.spec.js validates against exact data first.
//
// Each test imports /fieldmodels.js and /fieldlines2d.js in the page and runs a rule over a group
// of configurations in the page (the lines can be large, so only the list of problems comes back).
//
//   lectures  every `electrostatic_direct` example of tests/fixtures/examples.json, as an
//             electrostatic_line_charge model on the 668 x 335 canvas (normalised coordinates
//             x 668 / x 335; uniformX = uniformElecX / 335)
//   random    20 seeded random configurations of 3-4 charges of +-1..+-3
//   magnetic  the `magnetostatic_direct` examples, as magnetostatic_wire models
//
// These fail against the stubs in fieldlines2d.js and must pass once it is implemented. The tests
// use the plain Playwright `test` (not the fake-clock fixture of helpers.js), because the time
// budget test needs the real performance.now().

import { test, expect } from '@playwright/test';
import { examples } from './helpers.js';
import {
    lectureConfigs, randomConfigs, dragConfigs, rotationConfigs, axisConfigs,
} from './fieldline-configs.js';

const lectures = lectureConfigs(examples, 'electrostatic_direct');
const random = randomConfigs(20);
const magnetic = lectureConfigs(examples, 'magnetostatic_direct');
const randomMagnetic = randomConfigs(8, { seed: 4321 });
const fourCharge = randomConfigs(20).filter(c => c.charges.length === 4).slice(0, 5);

const MAX_REPORTED = 12;

/** Run rule `name` of fieldline-rules-lib.js in the page with the real modules; fail with a clear message. */
async function expectRule(page, name, args, title) {
    await page.goto('/tests/blank.html');
    const problems = await page.evaluate(async ({ name, args }) => {
        const [models, lines2d, rules] = await Promise.all([
            import('/fieldmodels.js'), import('/fieldlines2d.js'), import('/tests/fieldline-rules-lib.js'),
        ]);
        return rules[name]({ models, lines2d, now: () => performance.now() }, ...args);
    }, { name, args });
    if (problems.length) {
        const shown = problems.slice(0, MAX_REPORTED).map(p => `  - ${p}`).join('\n');
        const more = problems.length > MAX_REPORTED ? `\n  ... and ${problems.length - MAX_REPORTED} more` : '';
        throw new Error(`${title}: ${problems.length} violation(s)\n${shown}${more}`);
    }
    expect(problems).toEqual([]);
}

const electricGroups = [
    ['lectures', lectures],
    ['random configurations', random],
];

test.describe('electric field lines (2D line charges)', () => {
    for (const [group, configs] of electricGroups) {
        test.describe(group, () => {
            test('each source has exactly N|q| attached line ends (N lines per unit charge)', async ({ page }) => {
                await expectRule(page, 'ruleExactCounts', [configs], 'Exact counts');
            });
            test('every end is on a source (within 1 px) or at infinity (> 3 view diagonals away); from/to agree', async ({ page }) => {
                await expectRule(page, 'ruleEndpoints', [configs], 'Endpoints');
            });
            test('lines are tangent to the field (< 3 deg outside the drawn charges), oriented + to -', async ({ page }) => {
                await expectRule(page, 'ruleTangency', [configs], 'Tangency');
            });
            test('lines do not depend on the view: view x2 gives the same lines inside the original view (1 px)', async ({ page }) => {
                await expectRule(page, 'ruleViewportIndependence', [configs], 'Viewport independence');
            });
            test("Gauss's law: net signed crossings of random polygons = N x enclosed charge", async ({ page }) => {
                await expectRule(page, 'ruleGauss', [configs], "Gauss's law");
            });
            test('potential contours cross field lines at 90 +- 6 deg', async ({ page }) => {
                await expectRule(page, 'rulePerpendicularContours', [configs], 'Perpendicular contours');
            });
        });
    }

    test('potential contours do not depend on the view; contourStep is proportional to K', async ({ page }) => {
        await expectRule(page, 'ruleViewportIndependenceContours', [lectures, []], 'Potential contours: viewport independence');
        await expectRule(page, 'ruleContourStepScaling', [lectures], 'contourStep scaling');
    });

    test('dragging a charge in 1 px steps changes attachment angles by <= 3 deg per step and never the counts', async ({ page }) => {
        await expectRule(page, 'ruleDragStability', [dragConfigs()], 'Drag stability');
    });

    test('rotating a neutral configuration by 90 deg rotates the attachment angles (decision D5)', async ({ page }) => {
        await expectRule(page, 'ruleRotationInvariance', [rotationConfigs()], 'Rotation invariance');
    });

    test('straight lines along the axis of a pair exactly where a drawn level of psi lies (none for +-2 by default)', async ({ page }) => {
        await expectRule(page, 'ruleAxisLine', [axisConfigs()], 'Axis line');
    });
});

test.describe('magnetic field lines (straight wires)', () => {
    for (const [group, configs] of [['lectures', magnetic], ['random configurations', randomMagnetic]]) {
        test(`${group}: contours are closed without a uniform field, tangent to B, spaced by contourStep`, async ({ page }) => {
            await expectRule(page, 'ruleMagnetic', [configs], 'Magnetic lines');
        });
    }
    test('magnetic lines do not depend on the view', async ({ page }) => {
        await expectRule(page, 'ruleViewportIndependenceContours', [[], magnetic], 'Magnetic lines: viewport independence');
    });
});

test.describe('performance', () => {
    test('4-charge configurations: median of 10 calls < 30 ms (electric lines + potential contours; magnetic lines)', async ({ page }) => {
        await expectRule(page, 'ruleTimeBudget', [fourCharge], 'Time budget');
    });
});
