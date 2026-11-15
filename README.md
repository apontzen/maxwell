# MEME: an interactive electromagnetism visualiser

MEME draws electric and magnetic fields of charges and currents that you can drag around:
arrow ("quiver") plots, field lines, equipotentials, forces and torques, Gaussian surfaces.
It includes a simulation of Maxwell's equations in which moving charges radiate. It runs 
in the browser: the numerical core is Rust compiled to WebAssembly (`maxwell/`), and the 
drawing and user interface are plain JavaScript modules. 

It is licenced under BSD.

It was written in 2024 by Andrew Pontzen for Durham University [electromagnetism lectures](https://electromagnetism.pontzen.uk),
and you can access the latest published version here: [MEME](https://electromagnetism.pontzen.uk/MEME/).

While originally hand-coded, by 2026, Claude Code was capable of implementing this kind of physics 
in-browser, and the code was updated and improved. 

[HOW_IT_WORKS.md](HOW_IT_WORKS.md) explains the physics and the numerical methods (written
by Claude in 2026).

## Building

You need:

- Rust, with the WebAssembly target. Install it with `rustup`, then:

  ```bash
  rustup target add wasm32-unknown-unknown
  ```

- [`wasm-pack`](https://rustwasm.github.io/wasm-pack/):

  ```bash
  cargo install wasm-pack
  ```

- Node.js 18 or later (only for the development server and the tests).

Build the WebAssembly module into `maxwell/out/`:

```bash
npm run build:wasm
```

This is the same as `wasm-pack build --target web -d out` in `maxwell/`. Rebuild whenever you change
the Rust code. The JavaScript needs no build step.

The lecture notes are deployed with Rust 1.90 (see their workflow). Newer versions work too.

## Serving

The pages load ES modules and fetch the `.wasm` file, so they must be served over HTTP: opening
`index.html` directly from disk (`file://`) does not work. Any static file server pointed at the
repository root will do, for example the one used by the tests:

```bash
npm run serve
```

It serves the repository at <http://127.0.0.1:8123/> (pass another port as an argument:
`node tests/server.mjs 8000`). The pages are:

- `index.html`: the explorer, where you can add, edit and delete charges, choose the solver and
  copy or paste the state as JSON.
- `embed_test.html`: a page of embeds, as used in the lecture notes.

## Embedding

A page that includes the visualiser loads `ui.js` and calls `embed()`:

```html
<script type="module">
  import { embed } from './ui.js';
  embed();
</script>
```

Each `<meme-embed>` element on the page then becomes an interactive canvas. Its `meme` attribute holds
the starting state as JSON (the same format as the explorer's "Copy" button):

```html
<meme-embed fieldlines-checkbox="yes"
  meme='{"charges":[{"x":0.3,"y":0.5,"charge":2,"id":0},{"x":0.7,"y":0.5,"charge":-2,"id":1}],
         "solver":"electrostatic_direct","plotType":"fieldline"}'></meme-embed>
```

The state's keys:

| Key | Meaning |
|---|---|
| `charges` | List of `{x, y, charge, id}`, with `x` and `y` as fractions of the canvas width and height. Optional: `isTestCharge` (feels forces but makes no field), `dipoleWith` (the `id` of the partner in a rigid dipole). |
| `solver` | The field model, see below. |
| `plotType` | `"quiver"` (arrows) or `"fieldline"`. |
| `showPotential` | Draw equipotentials. |
| `perCharge` | Show each charge's field separately (quiver mode). |
| `uniformElecX`, `uniformElecY` | An external uniform field. |
| `dipoleMode` | Treat charges with `dipoleWith` as rigid dipoles. |
| `linesPerUnitCharge`, `fieldLineOffset` | Electric field-line density: lines per unit charge (default 5), and the offset of the drawn levels as a fraction of the line spacing (default 0.5). See HOW_IT_WORKS.md. |

Note that the `charges` are interpreted as currents directed through the +z axis (out of the screen) 
in the case of magnetic field solvers.

The solvers:

| `solver` | Model |
|---|---|
| `electrostatic_direct` | Point charges (3D field, shown in the plane). Field lines are drawn for line charges. |
| `electrostatic_line_charge` | Line charges perpendicular to the page (2D). |
| `electrostatic_fourier` | Line charges, with the field solved on a grid. |
| `dynamic` | Maxwell's equations on a grid: moving charges radiate. |
| `magnetostatic_moving_charge` | A snapshot of point charges moving through the page (3D). Field lines are drawn for straight wires. `magnetostatic_direct` is a deprecated alias. |
| `magnetostatic_wire` | Straight wires perpendicular to the page (2D). |

Optional attributes add controls to an embed: `fieldlines-checkbox`, `percharge-checkbox`,
`forces-checkbox`, `dipole-checkbox` (with a value such as `"yes"`), and `gaussian-surface`
(draw a closed curve with the mouse to see the flux through it).

## Testing

The regression tests (Playwright, in a real browser) and the Rust unit tests are described in
[tests/README.md](tests/README.md) and can be run with:

```bash
npm ci
```

```bash
npx playwright install chromium
```

```bash
npm test
```

```bash
cd maxwell && cargo test --release --lib
```
