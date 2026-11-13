// Shared helpers for the MEME regression tests.
//
// The tests deliberately interact with the visualiser only through the DOM, the
// canvas pixels and real mouse/touch events, never through internal JS state, so
// that they remain valid while the implementation is refactored.

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { test as base, expect } from '@playwright/test';

const here = path.dirname(fileURLToPath(import.meta.url));

/** All examples from the lecture notes and embed_test.html; see extract-examples.mjs */
export const examples = JSON.parse(fs.readFileSync(path.join(here, 'fixtures', 'examples.json'), 'utf8'));

/** A fixed instant at which the fake clock starts, so that runs are reproducible */
const CLOCK_START = new Date('2025-01-01T00:00:00Z');

/** How long a time-dependent ("dynamic") example keeps animating after an interaction */
export const DYNAMIC_RUNTIME_MS = 10_000;

/**
 * Test fixture that
 *  - installs a fake clock (so that requestAnimationFrame / performance.now are
 *    deterministic, which matters for the time-dependent solver), and
 *  - fails the test if the page throws or logs a console error.
 */
export const test = base.extend({
    page: async ({ page }, use, testInfo) => {
        const errors = [];
        page.on('pageerror', err => errors.push(`pageerror: ${err.message}`));
        page.on('console', msg => {
            if (msg.type() === 'error') errors.push(`console.error: ${msg.text()}`);
        });
        await installPausedClock(page);
        await use(page);
        if (!testInfo.annotations.some(a => a.type === 'allow-console-errors')) {
            expect(errors, 'page should not throw or log console errors').toEqual([]);
        }
    },
});
export { expect };

/**
 * Install a fake clock that only advances when the test says so (via
 * page.clock.runFor), making animation frames, timers and performance.now()
 * fully deterministic.
 */
export async function installPausedClock(page) {
    await page.clock.install({ time: CLOCK_START });
    await page.clock.pauseAt(new Date(CLOCK_START.getTime() + 1000));
}

// ---------------------------------------------------------------------------
// Example metadata
// ---------------------------------------------------------------------------

/** The checkboxes an embed should show, derived from its attributes and state. */
export function expectedCheckboxes(example) {
    const a = example.attributes;
    const s = example.state;
    const magnetic = s.solver === 'magnetostatic_direct';
    const boxes = [];
    // Order matches the order in which ui.js adds them to the DOM
    if (a['fieldlines-checkbox'])
        boxes.push({ key: 'fieldlines', label: 'Show field lines', checked: s.plotType === 'fieldline' });
    if (a['percharge-checkbox'])
        boxes.push({ key: 'percharge', label: magnetic ? 'Per current' : 'Per charge', checked: !!s.perCharge });
    if (a['forces-checkbox'])
        boxes.push({ key: 'forces', label: 'Show forces', checked: false });
    if (a['dipole-checkbox'])
        boxes.push({ key: 'dipole', label: magnetic ? 'Rigid current loop' : 'Rigid dipole', checked: !!s.dipoleMode });
    return boxes;
}

export function isDynamic(example) {
    return example.state.solver === 'dynamic';
}

export function hasGaussianSurface(example) {
    return example.attributes['gaussian-surface'] !== undefined;
}

/** Pairs [i, j] of charge indices that form rigid dipoles (each pair listed once). */
export function dipolePairs(example) {
    const charges = example.state.charges;
    const pairs = [];
    charges.forEach((c, i) => {
        if (c.dipoleWith !== undefined && c.dipoleWith > c.id) {
            const j = charges.findIndex(o => o.id === c.dipoleWith);
            if (j >= 0) pairs.push([i, j]);
        }
    });
    return pairs;
}

/**
 * Indices of charges that can be grabbed with the mouse, i.e. that are not hidden
 * underneath a later charge at the same position (later charges are drawn on top
 * and are picked first).
 */
export function grabbableCharges(example) {
    const charges = example.state.charges;
    return charges.map((_, i) => i).filter(i => !charges.slice(i + 1).some(
        o => Math.hypot(o.x - charges[i].x, (o.y - charges[i].y) / 2) < 0.03));
}

/**
 * Up to `max` grabbable charges, spread across the list, for examples with
 * many charges.
 */
export function sampleCharges(example, max = 4) {
    const idx = grabbableCharges(example);
    if (idx.length <= max) return idx;
    const out = new Set([idx[0], idx[Math.floor(idx.length / 2)], idx[idx.length - 1]]);
    return [...out].sort((a, b) => a - b);
}

/**
 * A deterministic drop target (normalised coordinates) for dragging a charge
 * from (x, y): move 0.2 of the canvas width/height towards the centre.
 */
export function dragTarget({ x, y }) {
    return { x: x < 0.5 ? x + 0.2 : x - 0.2, y: y < 0.5 ? y + 0.2 : y - 0.2 };
}

// ---------------------------------------------------------------------------
// Loading pages
// ---------------------------------------------------------------------------

export function harnessUrl(example) {
    return '/tests/harness.html?attrs=' + encodeURIComponent(JSON.stringify(example.attributes));
}

/** Load a single example in the harness page and wait until it has rendered. */
export async function gotoExample(page, example) {
    await page.goto(harnessUrl(example));
    await waitForRender(page, 1);
}

/**
 * Wait until `count` canvases matching `selector` exist and have been drawn on
 * (unless `allowBlank`), and the drawing is stable while the fake clock advances.
 */
export async function waitForRender(page, count, { selector = 'meme-embed canvas', allowBlank = false } = {}) {
    const deadline = Date.now() + 15_000;
    let previous = null;
    while (Date.now() < deadline) {
        // Advancing the fake clock lets timers and animation frames in the page run
        await page.clock.runFor(48);
        const hashes = await page.evaluate(({ selector, allowBlank }) => {
            const canvases = [...document.querySelectorAll(selector)];
            return canvases.map(c => {
                if (c.width === 0 || c.height === 0) return null;
                const data = c.getContext('2d').getImageData(0, 0, c.width, c.height).data;
                let nonBlank = 0, h = 0;
                for (let i = 0; i < data.length; i += 4) {
                    if (data[i + 3] !== 0) nonBlank++;
                    h = (h * 31 + data[i] + 7 * data[i + 1] + 13 * data[i + 2] + 17 * data[i + 3]) | 0;
                }
                return (nonBlank > 0 || allowBlank) ? h : null;
            });
        }, { selector, allowBlank });
        const ready = hashes.length === count && hashes.every(h => h !== null);
        if (ready && previous !== null && JSON.stringify(previous) === JSON.stringify(hashes)) {
            await alignToMidFrame(page);
            return;
        }
        previous = ready ? hashes : null;
        await page.waitForTimeout(20);
    }
    throw new Error(`Timed out waiting for ${count} canvas(es) matching "${selector}" to render`);
}

/**
 * Advance the fake clock to exactly half-way between two animation frames.
 *
 * The fake clock's origin, and hence the phase of the 16ms animation frames
 * relative to the start of a test, depends on page-load timing. The output of the
 * time-dependent solver depends on how user interactions interleave with frames,
 * so we fix the phase: run the clock until a frame fires, then move half a frame on.
 * Interactions then never coincide with frame boundaries.
 */
export async function alignToMidFrame(page) {
    await page.evaluate(() => {
        window.__memeTestFrameSeen = false;
        requestAnimationFrame(() => { window.__memeTestFrameSeen = true; });
    });
    for (let i = 0; i < 64; i++) {
        await page.clock.runFor(1);
        if (await page.evaluate(() => window.__memeTestFrameSeen)) {
            await page.clock.runFor(8);
            return;
        }
    }
    throw new Error('Fake animation frame never fired');
}

// ---------------------------------------------------------------------------
// Reading the canvas
// ---------------------------------------------------------------------------

export function embed(page, index = 0) {
    return page.locator('meme-embed').nth(index);
}

/**
 * The canvas of the index'th embed on the page. If `target` is already a Locator
 * (e.g. for the explorer page) it is returned unchanged.
 */
export function canvas(target, index = 0) {
    if (typeof target.nth === 'function') return target;
    return embed(target, index).locator('canvas');
}

/** The exact contents of the canvas backing store, as a PNG buffer. */
export async function canvasPng(page, index = 0) {
    const dataUrl = await canvas(page, index).evaluate(c => c.toDataURL('image/png'));
    return Buffer.from(dataUrl.split(',')[1], 'base64');
}

/** Raw RGBA pixel data of the canvas (for exact in-test comparisons). */
export async function canvasPixels(page, index = 0) {
    const { width, height, b64 } = await canvas(page, index).evaluate(c => {
        const bytes = new Uint8Array(c.getContext('2d').getImageData(0, 0, c.width, c.height).data.buffer);
        let s = '';
        for (let i = 0; i < bytes.length; i += 0x8000)
            s += String.fromCharCode.apply(null, bytes.subarray(i, i + 0x8000));
        return { width: c.width, height: c.height, b64: btoa(s) };
    });
    return { width, height, data: Buffer.from(b64, 'base64') };
}

/** Number of pixels that differ between two canvasPixels() results. */
export function countDifferentPixels(a, b) {
    if (a.width !== b.width || a.height !== b.height) return Infinity;
    let n = 0;
    for (let i = 0; i < a.data.length; i += 4) {
        if (a.data[i] !== b.data[i] || a.data[i + 1] !== b.data[i + 1] ||
            a.data[i + 2] !== b.data[i + 2] || a.data[i + 3] !== b.data[i + 3]) n++;
    }
    return n;
}

/** Average colour of the canvas in a small square around a normalised point. */
export async function colourAt(page, { x, y }, radius = 3, index = 0) {
    return canvas(page, index).evaluate((c, { x, y, radius }) => {
        const rect = c.getBoundingClientRect();
        const scale = c.width / c.clientWidth;
        const px = Math.round(x * rect.width * scale), py = Math.round(y * rect.height * scale);
        const r = Math.round(radius * scale);
        const d = c.getContext('2d').getImageData(px - r, py - r, 2 * r + 1, 2 * r + 1).data;
        const sum = [0, 0, 0, 0];
        for (let i = 0; i < d.length; i += 4) for (let k = 0; k < 4; k++) sum[k] += d[i + k];
        const n = d.length / 4;
        return sum.map(s => Math.round(s / n));
    }, { x, y, radius });
}

/** Rough classification of an RGBA colour, used for semantic checks. */
export function classifyColour([r, g, b, a]) {
    if (a < 128) return 'transparent';
    if (r > 180 && g < 80 && b < 80) return 'red';
    if (b > 180 && r < 80 && g < 80) return 'blue';
    return 'other';
}

// ---------------------------------------------------------------------------
// Interacting
// ---------------------------------------------------------------------------

/**
 * Convert normalised canvas coordinates (as used in the embed JSON) to page
 * coordinates. Charge positions are stored relative to the canvas'
 * getBoundingClientRect() size, and mouse offsets are relative to its padding box.
 */
export async function toPage(page, { x, y }, index = 0) {
    const m = await canvas(page, index).evaluate(c => {
        const r = c.getBoundingClientRect();
        return { left: r.left + c.clientLeft, top: r.top + c.clientTop, width: r.width, height: r.height };
    });
    // Round to whole pixels: mouse and touch events report coordinates with different
    // sub-pixel precision, which would otherwise make otherwise-identical drags differ.
    return { x: Math.round(m.left + x * m.width), y: Math.round(m.top + y * m.height) };
}

/**
 * Drag with the mouse between two normalised canvas positions. `frameMs` is the
 * amount of fake time that passes between successive mouse-move events (so that
 * animation frames run during the drag, as they would for a real user).
 */
export async function mouseDrag(page, from, to, { steps = 10, index = 0, frameMs = 0 } = {}) {
    const a = await toPage(page, from, index);
    const b = await toPage(page, to, index);
    await page.mouse.move(a.x, a.y);
    await page.mouse.down();
    for (let i = 1; i <= steps; i++) {
        await page.mouse.move(Math.round(a.x + ((b.x - a.x) * i) / steps), Math.round(a.y + ((b.y - a.y) * i) / steps));
        if (frameMs > 0) await page.clock.runFor(frameMs);
    }
    await page.mouse.up();
}

/** Draw a closed polygonal path (e.g. a Gaussian surface) with the mouse. */
export async function mouseDrawPath(page, points, { index = 0 } = {}) {
    const pts = [];
    for (const p of points) pts.push(await toPage(page, p, index));
    await page.mouse.move(pts[0].x, pts[0].y);
    await page.mouse.down();
    for (const p of pts.slice(1)) await page.mouse.move(p.x, p.y, { steps: 3 });
    await page.mouse.up();
}

/** Points on an ellipse in normalised coordinates (circle in CSS pixels for a 2:1 canvas). */
export function loopAround({ x, y }, radius = 0.08, n = 24, aspect = 2) {
    const pts = [];
    for (let i = 0; i <= n; i++) {
        const t = (2 * Math.PI * i) / n;
        pts.push({ x: x + radius * Math.cos(t), y: y + radius * aspect * Math.sin(t) });
    }
    return pts;
}

/** Drag with a single finger, using real touch events via the Chrome DevTools Protocol. */
export async function touchDrag(page, from, to, { steps = 10, index = 0 } = {}) {
    const a = await toPage(page, from, index);
    const b = await toPage(page, to, index);
    const cdp = await page.context().newCDPSession(page);
    const point = (x, y) => [{ x, y, radiusX: 1, radiusY: 1, force: 1, id: 0 }];
    await cdp.send('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: point(a.x, a.y) });
    for (let i = 1; i <= steps; i++) {
        const x = Math.round(a.x + ((b.x - a.x) * i) / steps);
        const y = Math.round(a.y + ((b.y - a.y) * i) / steps);
        await cdp.send('Input.dispatchTouchEvent', { type: 'touchMove', touchPoints: point(x, y) });
    }
    await cdp.send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
    await cdp.detach();
}

/** Let the page process pending work (and animation frames) for `ms` of fake time. */
export async function advance(page, ms = 100) {
    await page.clock.runFor(ms);
}

export function checkbox(page, key, index = 0) {
    return embed(page, index).locator(`.${key}-controls input[type=checkbox]`);
}
