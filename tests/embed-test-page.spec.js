// Tests for embed_test.html as a whole page: several embeds on one page share
// the WebAssembly module and must initialise and respond independently.

import {
    test, expect, examples, expectedCheckboxes, waitForRender, embed, canvasPng, canvasPixels,
    countDifferentPixels, mouseDrag, dragTarget, advance, checkbox,
} from './helpers.js';

const pageExamples = examples.filter(e => e.source.startsWith('embed_test.html'));

test.describe('embed_test.html', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/embed_test.html');
        await waitForRender(page, pageExamples.length);
    });

    test('all embeds render', async ({ page }) => {
        await expect(page.locator('meme-embed')).toHaveCount(pageExamples.length);
        for (const [i, example] of pageExamples.entries()) {
            expect(await canvasPng(page, i)).toMatchSnapshot([`embed-${i + 1}.png`]);
            const expected = expectedCheckboxes(example);
            await expect(embed(page, i).locator('input[type=checkbox]')).toHaveCount(expected.length);
            for (const cb of expected) {
                await expect(embed(page, i).locator(`.${cb.key}-controls label`)).toHaveText(cb.label);
                if (cb.checked) await expect(checkbox(page, cb.key, i)).toBeChecked();
                else await expect(checkbox(page, cb.key, i)).not.toBeChecked();
            }
        }
    });

    test('interacting with one embed leaves the others untouched', async ({ page }) => {
        const before = [];
        for (let i = 0; i < pageExamples.length; i++) before.push(await canvasPixels(page, i));

        // drag a charge in the first embed and toggle the tickbox in the second
        const c = pageExamples[0].state.charges[0];
        await page.locator('meme-embed').nth(0).scrollIntoViewIfNeeded();
        await mouseDrag(page, c, dragTarget(c), { index: 0 });
        await checkbox(page, 'fieldlines', 1).click();
        await advance(page);

        expect(await canvasPng(page, 0)).toMatchSnapshot(['embed-1-dragged.png']);
        expect(await canvasPng(page, 1)).toMatchSnapshot(['embed-2-fieldlines.png']);
        for (let i = 2; i < pageExamples.length; i++) {
            expect(countDifferentPixels(before[i], await canvasPixels(page, i)), `embed ${i + 1} changed`).toBe(0);
        }
    });
});
