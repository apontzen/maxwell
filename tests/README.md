# Regression tests

Browser tests (Playwright + Chromium) intended as a safety net while refactoring.
They pin down the current behaviour of every interactive example used in the
PHYS1122 lecture notes (`apontzen/du_electromagnetism`), of `embed_test.html`, and
of the stand-alone explorer (`index.html`).

The tests drive the visualiser only through the DOM, real mouse/touch events
and the canvas pixels. They never touch internal JS state, so they stay valid
while the internals change.

## Running

Prerequisites: Rust with the `wasm32-unknown-unknown` target, `wasm-pack`, Node ≥ 18.

```bash
npm run build:wasm
```

```bash
npm ci
```

```bash
npx playwright install chromium
```

```bash
npm test
```

Useful variants:

```bash
npx playwright test -g lecture11-2
```

```bash
npx playwright test tests/explorer.spec.js
```

```bash
npx playwright test --ui
```

The test runner starts its own static server (`tests/server.mjs`, port 8123).
The full suite takes about 2 minutes on 4 cores.

### Reference images and platforms

The reference images in `tests/__snapshots__/` are for **Linux** (they were
generated with Playwright 1.63.0 / Chromium 1243, as used in CI). Canvas
anti-aliasing can differ on macOS, so on a Mac run the suite in the same
container as CI:

```bash
docker run --rm --ipc=host -v "$PWD":/work -w /work mcr.microsoft.com/playwright:v1.63.0-noble /bin/bash -c "npm ci && npx playwright test"
```

(`maxwell/out` must already have been built on the host; the WebAssembly is
platform-independent.) Alternatively, check out the pre-refactor commit, run
`npm run test:update` to create local `*-darwin.png` baselines (these are
git-ignored), then switch back to the refactor branch and run `npm test`.

## When something fails

Failures save the actual image and a diff under `test-results/`. If a change is
*intended* (e.g. a deliberate restyle), regenerate the references and review the
image changes in the git diff before committing:

```bash
npm run test:update
```

## What is tested

`embeds.spec.js` runs the following for each example in `fixtures/examples.json`
(23 from the lecture notes and 7 from `embed_test.html`). Each example is loaded
on its own in `harness.html`, which copies the lecture site's CSS.

| Test | Checks |
|---|---|
| initial render | Reference image of the canvas; canvas size; exactly the expected tickboxes, with the right labels ("Rigid dipole" vs "Rigid current loop", …) and initial states; each electric charge drawn in the right colour at its position |
| toggle tickbox | Reference image after toggling; the image changes; toggling back restores the original exactly |
| drag charge *i* | Reference image after a mouse drag (up to 4 charges per example); the charge appears at the drop point and has left its original position |
| dynamic solver drags | Radiation from a dragged charge at t = 0, 0.5 s and 2 s; the animation stops after 10 s |
| touch drag | A touch drag gives exactly the same result as the mouse drag |
| rigid dipoles | Rotating by dragging an end, translating by dragging the midpoint, and free movement with the tickbox off |
| Gaussian surfaces | Drawing a surface around one charge, around all charges, around nothing; replacing it; dragging a charge into it |
| click / double-click | Embeds don't allow editing, so clicks change nothing and no properties panel appears |
| redraw | Redrawing an unchanged state reproduces the first render exactly |
| phone width / resize | Rendering at 390 px wide, and after resizing the window down and back up |

`embed-test-page.spec.js` loads `embed_test.html` with its seven embeds and
checks that they all initialise (they share one WebAssembly instance) and that
interacting with one leaves the others untouched.

`explorer.spec.js` covers `index.html`: adding charges (button and double-click),
the properties panel (strength slider, test-charge tickbox, delete), dragging,
clear, every solver with its control visibility, the potential, field-line,
force and per-charge tickboxes under both static solvers, the dipole tickbox, rigid dipoles, persistence in
localStorage across reloads, and copy/paste of the JSON state.

`fixtures.spec.js` checks that `fixtures/examples.json` still matches the
lecture notes. It runs whenever a `du_electromagnetism` checkout is found:
alongside this repo, as its parent (when this repo is its `MEME` submodule), or
at `$DU_ELECTROMAGNETISM`. It is skipped otherwise. When the notes change:

```bash
npm run extract-examples
```

```bash
npm run test:update
```

## Determinism

* Playwright's fake clock is installed and paused, so `requestAnimationFrame`,
  `performance.now()` and timers only advance when a test says so
  (`page.clock.runFor`). This makes the time-dependent solver reproducible.
* The fake clock's origin depends on page-load timing, and the dynamic solver's
  output depends on how interactions interleave with animation frames. So after
  each render-wait the clock is moved to exactly half-way between two frames
  (`alignToMidFrame` in `helpers.js`).
* Images come from the canvas backing store (`toDataURL`), not from screenshots,
  so overlaid controls and fonts don't affect them. `deviceScaleFactor` is 1.
* Pointer coordinates are rounded to whole pixels, because mouse and touch
  events report sub-pixel positions differently.
* WebAssembly built with Rust 1.90 (as used to deploy the lecture notes) or 1.97,
  with or without `wasm-opt`, passes against the same reference images.

## Known bugs (tests marked `test.fail`)

These tests are *expected* to fail on the current code. When a bug is fixed,
Playwright reports the test as an unexpected pass. Then remove its `test.fail`
(and, for the second bug, the `redrawDiffersFromFirstRender` special case in
`helpers.js`).

1. **Double-clicking a charge throws `ReferenceError: deleteCharge is not defined`**
   in both embeds and the explorer. `deleteCharge` is declared inside an `if`
   block in `ui.js`, and in (strict) module code that makes it block-scoped. In
   the explorer this means double-click-to-delete does nothing.
2. **Field-line direction arrows inherit the previous frame's colour** in
   magnetostatic field-line plots. `drawDirectionArrow` in `fieldline.js` never
   sets `strokeStyle`. With rigid current loops shown, the arrowheads are black
   on first render and turn grey or purple after any interaction.
3. **"Charge"/"Current" labels never switch in the explorer.** `ui.js` queries
   `.charge-or-current`, but `index.html` uses the class `charge_or_current`.

A further quirk is recorded as-is rather than marked as a bug: after resizing
the window down and back up, the canvas comes back 1 px taller than it started.
