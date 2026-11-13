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
force and per-charge tickboxes, the dipole tickbox, rigid dipoles, persistence in
localStorage across reloads, and copy/paste of the JSON state. It also checks
that ticking field lines visibly switches a 3D solver to its 2D counterpart and
back, and that the deprecated `magnetostatic_direct` name renders exactly as
`magnetostatic_moving_charge`.

### Field lines (2D)

Field lines are always drawn for a 2D world (line charges / straight wires, or the 2D grid
solvers); see [issue #3](https://github.com/apontzen/maxwell/issues/3). Their geometry
(`fieldlines2d.js`) is tested directly, independently of the pictures:

- `fieldline-checks.js` contains pure geometric checkers.
- `fieldline-checks.spec.js` unit-tests those checkers against exact analytic
  field lines (for a ±q pair of line charges, the field lines are circular arcs).
  It also checks that deliberately broken implementations are caught.
- `fieldlines-rules.spec.js` applies the rules to every lecture example and to
  seeded random 3–4 charge configurations:
  - exactly 4|q| line ends per charge;
  - lines end only on charges or at infinity, and are tangent to the field;
  - lines and contours are independent of the view;
  - lines move continuously under 1 px drags;
  - rotation invariance and a line along the axis (decision D5);
  - equipotentials are perpendicular to field lines;
  - net crossings of any closed loop equal 4 × the enclosed charge;
  - magnetic field lines are closed;
  - a time budget per frame.
- `fieldlines-grid.spec.js` does the same for the grid solvers (`electrostatic_fourier`,
  `dynamic`; `gridfield.js` on top of the stream function computed in `maxwell/src/streamfn.rs`):
  4|q| line ends per charge, lines end on a charge, on the edge of the canvas (the grid world ends
  there) or fade out where there is no field, tangency, Gauss's law, a time budget, and closed
  field lines radiated by a shaken dipole. The lecture examples, random configurations and a net
  charge are used for the Fourier solver; for the dynamic solver, lecture 14 before and after a
  drag and lecture 16 after shaking a charge.
- The Rust unit tests (`cargo test --lib` in `maxwell/`) check the stream function's least-squares
  solver, and that it reproduces the grid's face fluxes exactly for pure radiation, closely
  (outside the smeared charges) for a static field, and continuously as a charge crosses cells.
  They also check the dynamic solver: static charges keep their field exactly, Gauss's law holds
  during sub-cell drags, the field settles to the static one after a drag, and radiation leaves
  through the absorbing layer.
- `solvers.spec.js` checks the solver registry (`solvers.js`): which model is used in
  each mode, the deprecated alias, and that forces and fields come from the same 2D model.

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

## Known quirks

One quirk is recorded as-is rather than treated as a bug: after resizing
the window down and back up, the canvas comes back 1 px taller than it started.
