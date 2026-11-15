# How MEME works

This note explains the physics and the mathematics behind what MEME draws: the field models, how
field lines are constructed, and how the time-dependent simulation of Maxwell's equations works. It
assumes vector calculus, Maxwell's equations in differential form and some familiarity with
numerical methods (finite differences, Fourier transforms), at the level of a second- or third-year
physics undergraduate. File names in brackets point to where each idea is implemented.

## 1. Conventions and units

Positions are in canvas pixels, with $x$ to the right and $y$ **downwards**, as on a screen. A
right-handed coordinate system then has $\hat z$ pointing *into* the screen. That fixes the sign of
$B_z$ and of the stream function below; nothing else depends on it.

Units are arbitrary: each model has its own overall constant, chosen so that the pictures have
sensible sizes. Field-line pictures don't depend on overall constants at all.

Fields shown as arrows are *softened* near each charge: $r^2 \to r^2 + s$ with $s = 5\,\text{px}^2$,
so that they stay finite at the charge. Field lines use the exact, unsoftened fields.

## 2. The field models

| Solver key | Sources | Field |
|---|---|---|
| `electrostatic_direct` | point charges in the plane of the screen | $\vec E = F \sum_i q_i\, \vec r_i / r_i^3$ (3D Coulomb, evaluated in the plane) |
| `electrostatic_line_charge` | infinite line charges perpendicular to the screen | $\vec E = K \sum_i q_i\, \vec r_i / r_i^2$ (2D Coulomb) |
| `magnetostatic_moving_charge` | point charges crossing the screen (a snapshot) | 3D Biot–Savart: $\vec B = F \sum_i q_i\, \vec r_i^{\,\perp} / r_i^3$ |
| `magnetostatic_wire` | infinite straight wires perpendicular to the screen | $\vec B = K \sum_i I_i\, \vec r_i^{\,\perp} / r_i^2$, i.e. $\lvert B\rvert \propto I/R$ |
| `electrostatic_fourier` | line charges, smeared over a few grid cells | 2D Poisson equation solved on a grid (section 4) |
| `dynamic` | moving line charges | 2D Maxwell's equations on a grid (section 5) |

Here $\vec r_i = \vec r - \vec r_i^{\,\text{src}}$, $r_i = |\vec r_i|$, and
$\vec r^{\,\perp} = (y, -x)$ is $\vec r$ rotated by $90^\circ$: magnetic field lines circle the
sources. Any of the models can include a uniform external field.

**Normalisation.** The 3D constant is $F = 20000$. The 2D one is chosen so that a line charge and a
point charge give the same field strength at a reference distance $r_\text{ref}$ of a quarter of the
canvas width $W$:

$$K = F / r_\text{ref}, \qquad r_\text{ref} = W/4.$$

Uniform fields are the same in both. A magnetic model is the corresponding electric one rotated by
$90^\circ$ in the plane: a wire carrying current $I$ has the field of a line charge $I$, rotated.

### 2.1 Why field lines are always drawn in 2D

In field-line mode, MEME draws the field of line charges (or wires) even when point charges
(or moving charges) are chosen. The arrows, forces and everything else in that mode switch to the
same 2D model, so the picture stays consistent.

The reason is Gauss's law. We would like a field-line picture to satisfy three rules:
- each charge $q$ has a number of lines proportional to $q$;
- lines start and end only on charges (or at infinity);
- the density of lines is proportional to $\lvert \vec E\rvert$.

In a plane these can all hold only if the field has zero divergence in that plane away from the
charges: the number of lines crossing a closed curve is then the flux through it, which is
proportional to the enclosed charge. For a line charge, $\nabla_{2\text{D}} \cdot (\vec r / r^2) = 0$
away from the charge, which is just 2D Gauss's law. But the field of a *point* charge, restricted to
a plane through it, has

$$\nabla_{2\text{D}} \cdot \frac{\vec r}{r^3} = \frac{2}{r^3} - \frac{3}{r^3} = -\frac{1}{r^3} \neq 0,$$

because flux leaks out of the plane. No set of lines in the plane can then satisfy all three rules
at once: if the lines are evenly spaced around each charge, their density falls off as $1/r$, not
$1/r^2$. (This is the point of the footnote in Lecture 1 of the notes.)

## 3. Field lines in 2D

### 3.1 The stream function

In 2D, a field with zero divergence (away from sources) can be written in terms of a scalar
*stream function* $\psi$. For the electric field of line charges, MEME uses

$$\vec E = 2\pi K \left( \frac{\partial \psi}{\partial y},\; -\frac{\partial \psi}{\partial x} \right),
\qquad
\psi = \frac{1}{2\pi} \sum_i q_i\, \theta_i + \frac{u_x (y - y_c) - u_y (x - x_c)}{2\pi K},$$

where $\theta_i$ is the polar angle of $\vec r$ about charge $i$, $\vec u$ is the uniform field and
$(x_c, y_c)$ is the centre of the canvas. Two properties make $\psi$ the right object:

1. **$\vec E$ is tangent to the contours of $\psi$**, because $\vec E \cdot \nabla \psi = 0$. So
   field lines are contours $\psi = \text{const}$.
2. **$\psi$ measures flux.** The flux of $\vec E$ across any curve from $\vec a$ to $\vec b$ is

   $$\int_{\vec a}^{\vec b} \vec E \cdot \hat n \, \mathrm{d}l = 2\pi K \left[ \psi(\vec b) - \psi(\vec a) \right]$$

   (with $\hat n$ on one fixed side of the direction of travel), independent of the path between
   them, provided the path doesn't wind round a charge.

Going once round a charge $q$, $\theta$ increases by $2\pi$, so $\psi$ increases by $q$: the flux
out of the charge is $2\pi K q$, as Gauss's law requires. This also means $\psi$ is **multivalued**.
It is defined only up to multiples of the charges, and computing it along a path means following
each $\theta_i$ continuously ("unwrapping") rather than taking its principal value.

### 3.2 Choosing which lines to draw

MEME draws the contours

$$\psi = \frac{k + f}{N}, \qquad k \in \mathbb{Z},$$

with $N$ lines per unit charge (default 5) and an offset $f$ (default $\tfrac12$). These are
configurable per embed.

Because adjacent contours differ by $\Delta\psi = 1/N$, each strip between adjacent lines carries
the same flux, $2\pi K / N$. Several things follow exactly:

- **Line count.** A charge $q$ has exactly $N\lvert q\rvert$ lines (if that is a whole number).
- **Line density.** The spacing of the lines is $\Delta\psi / \lvert\nabla\psi\rvert = 2\pi K / (N \lvert\vec E\rvert)$,
  so their density is exactly proportional to $\lvert \vec E\rvert$.
- **Gauss's law, by counting lines.** Going round any closed curve, $\psi$ increases by the
  enclosed charge $Q$. The number of levels it passes upwards, minus the number it passes downwards,
  is therefore $NQ$: lines crossing outwards minus lines crossing inwards equals $N$ times the
  enclosed charge. The Gaussian-surface tool in Lecture 5 relies on this.

**The offset $f$.** The constant in $\psi$ is a free choice, but it decides where the lines go.
- **Along a dipole's axis:** on the line through a pair of charges, outside both, $\psi = 0$. With
  $f = 0$ that is always a drawn level, so a dipole had a straight line running out along its axis
  in each direction. With $f = \tfrac12$ it isn't.
- **Translations:** moving every charge by the same amount moves the lines with them, whatever $f$
  is, since $\psi$ depends only on the angles $\theta_i$.
- **Rotations:** rotating the whole configuration by $\alpha$ changes every $\theta_i$ by $\alpha$,
  so $\psi$ changes by $Q_\text{total}\, \alpha / 2\pi$. For a neutral configuration the lines simply
  rotate with it.

### 3.3 Tracing a contour

A contour of $\psi$ is traced from a starting point by repeated steps
(`levelsets2d.js`, `traceContour`):

1. **Predictor.** Take a step of length $h$ along the unit tangent
   $\hat t = (\psi_y, -\psi_x)/\lvert\nabla\psi\rvert$, evaluated at the midpoint of the step.
2. **Corrector.** Project back onto the contour with Newton's method in the direction of the
   gradient:

   $$\vec x \leftarrow \vec x - \frac{\psi(\vec x) - L}{\lvert\nabla\psi\rvert^2}\, \nabla\psi,$$

   which converges quadratically to the level $L$.

Because of the corrector, errors don't accumulate: a long line stays exactly on its contour (to
rounding error), so it ends exactly on the right charge, and two lines never drift into each other.
The step length adapts to the curvature: a straight chord of length $h$ across an arc of curvature
$\kappa$ deviates from it by about $\kappa h^2/8$, which is kept below 0.1 px where the line is
visible.

Lines that leave the view are followed until they end on a charge or reach "infinity" (ten view
diagonals away). So whether two pieces of line seen on screen are connected doesn't depend on the
size of the view.

### 3.4 Where lines start

Every contour must be found exactly once.

- **Lines from charges.** Close to a charge $i$, $\psi \approx q_i \theta_i / 2\pi + \text{const}$ is
  monotonic in the angle round it. So on a small circle round the charge, $\psi$ passes through each
  level once. The $N\lvert q_i\rvert$ crossings are found by bracketing and a root-finder; each is
  the start of a line.
  - Lines are traced from positive charges along $\vec E$ until they reach a negative charge,
    where they use up one of its crossings.
  - Crossings of negative charges left unused belong to lines coming from infinity, and are traced
    backwards.
- **Lines that touch no charge.** With a uniform field, some lines touch no charge at all. Each of
  these crosses the boundary of the view, so they are found from the levels' crossings of the
  boundary. Seeds on a line already traced are skipped.
- **Closed loops.** Possible only in the time-dependent case (section 6.3).

### 3.5 Stagnation points

Where $\vec E = 0$, the contours of $\psi$ meet at a saddle. A contour exactly at a saddle's level
is an "X", and a tracer arriving there must decide which way to go.

- **Finding them.** For line charges, the stagnation points can be found exactly. With
  $z = x + iy$, the conjugate field is an analytic function,

  $$E_x - i E_y \propto \sum_i \frac{q_i}{z - z_i} + \frac{\bar u}{K},$$

  so its zeros are the roots of a polynomial of degree $n$ (for $n$ charges).
- **Passing through.** When a contour runs into a saddle, the tracer leaves along the next branch
  in a fixed rotational order. This pairs the four branches consistently, so each is drawn once.

### 3.6 Magnetic field lines and equipotentials

For wires, $\vec B = \nabla \times (A_z \hat z) = (\partial_y A_z, -\partial_x A_z)$, with the
single-valued vector potential

$$A_z = K \sum_i I_i \ln r_i + (\text{uniform field term}).$$

Magnetic field lines are contours of $A_z$, spaced by a fixed step $\Delta A$. The flux of $\vec B$
between adjacent lines (per unit length of wire) is $\Delta A$, so again density $\propto \lvert
\vec B\rvert$.

Because $A_z$ is single-valued (there are no magnetic charges), every line closes: the magnetic
Gauss's law $\nabla\cdot\vec B = 0$ in pictures. Near a wire, $A_z \approx K I \ln r$, so a step
$\Delta A = \pi K/4$ gives $4\lvert I\rvert$ loops for every factor $e^\pi \approx 23$ in radius.
A wire can't have a fixed total number of loops: the flux diverges logarithmically both at the wire
and far away.

Equipotentials of line charges are contours of the 2D potential
$\phi = -K \sum_i q_i \ln(r_i / W) - \vec u \cdot (\vec r - \vec r_c)$, with the same step. They
are drawn red where $\phi > 0$ and blue where $\phi < 0$, and cross the field lines at right
angles, since $\nabla\phi \propto -\vec E \perp \nabla\psi$.

## 4. Grid solvers: the static field

The `electrostatic_fourier` and `dynamic` solvers work on a grid. This lets them handle situations
with no closed-form solution: the dynamic one needs a grid to evolve the fields at all.

### 4.1 The grid

The grid covers the canvas with square cells of side $\Delta$, about 6 px. It is surrounded by an
absorbing layer 20 cells thick (section 5.3). Quantities live at different points of each cell, a
staggered (Yee) arrangement:

| Quantity | Position in cell $(i, j)$, which spans $[i, i+1] \times [j, j+1]$ |
|---|---|
| $\rho$, $\phi$ | centre, $(i+\tfrac12, j+\tfrac12)$ |
| $E_x$ | middle of the right face, $(i+1, j+\tfrac12)$ |
| $E_y$ | middle of the bottom face, $(i+\tfrac12, j+1)$ |
| $B_z$ | corner, $(i+1, j+1)$ |

With this arrangement the natural finite differences are centred, and the discrete divergence of
a discrete curl vanishes identically. The discrete Gauss's law then behaves exactly like the
continuous one (section 5.2).

Each charge is spread over the grid as a Gaussian of width $\sigma = 2$ cells, centred on its exact
position. The grid's field then changes smoothly as a charge is dragged by less than a cell, and the
field of a single grid cell (where differencing errors are largest) is never shown directly.

### 4.2 The static field in free space

The static field is $\vec E = -\nabla\phi$ with the 5-point discrete Poisson equation

$$\frac{\phi_{i+1,j} - 2\phi_{i,j} + \phi_{i-1,j}}{\Delta^2} + \frac{\phi_{i,j+1} - 2\phi_{i,j} + \phi_{i,j-1}}{\Delta^2} = -\rho_{i,j}$$

(up to the code's sign conventions). A Poisson equation needs boundary conditions. The natural
choice for a grid, periodic boundaries solved by FFT, puts the charges in an infinite periodic array
of copies of themselves. The copies of a dipole distort its field visibly in a canvas only twice as
high as the dipole is long.

Instead MEME imposes **vacuum** boundary conditions (`freespace.rs`):

1. **The boundary ring.** The outermost ring of cells is at least 20 cells from any charge. There
   the potential of the smeared charges equals that of point charges,
   $\phi = -\sum_i Q_i \ln r_i / 2\pi$, to within $\sim e^{-r^2/2\sigma^2}$ (negligible). The
   potential on the ring is set to that.
2. **The interior.** It is then solved *exactly* for the discrete equation with these Dirichlet
   values. The eigenvectors of the discrete Laplacian with Dirichlet boundaries are sines, so a
   discrete sine transform in each direction diagonalises it, with eigenvalues
   $-\tfrac{4}{\Delta^2}\left[\sin^2\tfrac{\pi k}{2(n_x+1)} + \sin^2\tfrac{\pi l}{2(n_y+1)}\right]$.

The result agrees with the field of line charges in vacuum to better than $0.5^\circ$ in direction,
beyond 40 px from the charges. It also satisfies the grid's Gauss's law exactly, which the dynamic
solver needs: it starts from this field.

## 5. The dynamic solver: Maxwell's equations on a grid

### 5.1 The equations

With line charges and currents in the plane, the fields $E_x$, $E_y$ (in the plane) and $B_z$
(perpendicular to it) form a closed system, the "transverse electric" (TE) polarisation. In units
with $\varepsilon_0 = \mu_0 = c = 1$:

$$\frac{\partial B_z}{\partial t} = \frac{\partial E_x}{\partial y} - \frac{\partial E_y}{\partial x},
\qquad
\frac{\partial E_x}{\partial t} = \frac{\partial B_z}{\partial y} - J_x,
\qquad
\frac{\partial E_y}{\partial t} = -\frac{\partial B_z}{\partial x} - J_y,$$

with the constraint $\nabla\cdot\vec E = \rho$ (Gauss) and $\nabla\cdot\vec B = 0$ automatically
satisfied.

- **Faraday's law** gives the first equation; **Ampère–Maxwell** gives the other two.
- **The speed of light** is 1 px per unit of simulation time. One unit is 5 ms of real time, so
  light travels 200 px per second on screen, slow enough to watch.

These are discretised on the staggered grid with the **leapfrog** (Yee) scheme:
- $B_z$ lives half a time step behind $\vec E$;
- each update uses centred differences in space and time, so it is second-order accurate and
  conserves a discrete energy.

The scheme is stable if $c\,\Delta t \le \left(\Delta x^{-2} + \Delta y^{-2}\right)^{-1/2}$ (the
Courant condition). MEME takes 0.3 of that limit, because larger steps make the absorbing layer
reflect more.

### 5.2 Moving charges and Gauss's law

Gauss's law is not evolved; it is a constraint. Take the divergence of the Ampère–Maxwell equations:

$$\frac{\partial}{\partial t}\left(\nabla\cdot\vec E - \rho\right) = -\left(\nabla\cdot\vec J + \frac{\partial\rho}{\partial t}\right).$$

So if Gauss's law holds initially, it keeps holding exactly when the current satisfies the
continuity equation. On the staggered grid the same identity holds for the discrete operators.

So when a charge is dragged, MEME constructs the current density $\vec J$ from the change in its
smeared charge distribution, such that the *discrete* continuity equation holds exactly
(Esirkepov's method, `make_currents`). The charge is a product of 1D weights $W_x W_y$; its change
over a step is split into an $x$ part and a $y$ part,

$$\Delta\rho = q\left[\Delta W_x\, \frac{W_y^\text{old} + W_y^\text{new}}{2} + \frac{W_x^\text{old} + W_x^\text{new}}{2}\, \Delta W_y\right],$$

and each part is the discrete divergence of $J_x$ or $J_y$, which is then a running sum along that
direction.

With this current, the grid's Gauss's law is preserved to rounding error for as long as the
simulation runs, however the charges move. The flux out of a charge always equals the charge.

### 5.3 The absorbing boundary

Radiation reaching the edge of the grid must leave without reflecting, as if space continued.
MEME uses a **perfectly matched layer** (PML, Berenger's split-field form; `absorber.rs`).

**The PML.** In the layer, $B_z$ is split into two parts that are damped separately:

$$\frac{\partial B_{zx}}{\partial t} = -\frac{\partial E_y}{\partial x} - \sigma_x B_{zx},
\quad
\frac{\partial B_{zy}}{\partial t} = \frac{\partial E_x}{\partial y} - \sigma_y B_{zy},
\quad
\frac{\partial E_x}{\partial t} = \frac{\partial B_z}{\partial y} - \sigma_y E_x,
\quad
\frac{\partial E_y}{\partial t} = -\frac{\partial B_z}{\partial x} - \sigma_x E_y,$$

with $B_z = B_{zx} + B_{zy}$.

- **Why it doesn't reflect.** At a single frequency $\omega$, these equations are Maxwell's
  equations with $x$ replaced by a complex coordinate,
  $\partial_x \to (1 + i\sigma_x/\omega)^{-1}\partial_x$. A wave crossing into the layer keeps its
  form but decays as $\exp\!\left(-\cos\theta \int \sigma_x \,\mathrm{d}x\right)$. Its
  impedance is unchanged, so in the continuum there is no reflection at the interface, at any angle
  of incidence and any frequency.
- **On the grid.** A sudden jump in $\sigma$ would reflect, so $\sigma$ rises gradually, as
  $(\text{depth}/\text{thickness})^3$. Its maximum is chosen so that a wave entering the layer at
  normal incidence, reflecting off the outer edge and coming back is reduced in amplitude by a
  factor $10^6$. The damping is integrated semi-implicitly (Crank–Nicolson), which is
  unconditionally stable.

**Keeping the static field.** A PML on its own would also drain away the static field of the
charges, whose field lines reach into the layer. So the layer acts only on the *deviation* of the
fields from a static target field $\vec T$: $\vec E - \vec T$ is damped, not $\vec E$.

**Causality of the target.** The target is the static field of each charge where it was when the
light now arriving at that cell left it: the gradient of a retarded potential,

$$\phi_\text{ret}(\vec r, t) = \sum_i q_i\, G\big(\vec r - \vec r_i(t_\text{ret})\big),
\qquad t_\text{ret} = t - \lvert\vec r - \vec r_i\rvert / c,$$

where $G$ is the free-space potential of a unit smeared charge on the grid. The target is also
low-pass filtered (time constant 1 s), so that the layer enforces only the slowly varying part of
the field, not the radiation, which it must absorb. So the layer, like everything else, only learns
that a charge has moved when light from it arrives.

**The steady state.** When nothing has moved for longer than light takes to cross the grid, the
target is exactly the free-space static field. That is then the steady state of the whole system:
the field relaxes to that of charges in vacuum. In tests, after a drag the remaining difference
from a simulation on a grid five times larger is below 1% of the transient.

## 6. Field lines of a grid field

Field lines of the grid solvers (`streamfn.rs`, `gridfield.js`) are again contours of a stream
function, but $\psi$ must now be computed from the grid's field.

### 6.1 Splitting off the charges

Near a charge, $\psi$ winds by $q$ and is singular, which a grid can't represent. So

$$\psi = \psi_\text{sing} + \psi_\text{rest}, \qquad \psi_\text{sing} = \frac{1}{2\pi}\sum_i q_i\,\theta_i,$$

with $\psi_\text{sing}$ evaluated analytically about each charge's exact position. The remainder
$\psi_\text{rest}$ is single-valued and smooth, and is stored on the grid's corners.

It is found from the grid's own fluxes. The flux $E_x \Delta y$ through each face of the grid is
known exactly; the singular part's flux through the same face is subtracted, leaving residuals
$r_f$. $\psi_\text{rest}$ is then the least-squares solution of

$$\min_{\psi_\text{rest}} \sum_{\text{faces } f} \big(\delta_f \psi_\text{rest} - r_f\big)^2,$$

where $\delta_f$ is the difference of $\psi_\text{rest}$ across face $f$. Its normal equations are
a discrete Poisson equation with Neumann boundaries, diagonalised by a discrete cosine transform.

- **No charge present:** the grid field is exactly divergence-free, so this reproduces every face
  flux exactly. This is the case for pure radiation, e.g. after a shaken dipole has come back to
  rest on top of itself.
- **Static field:** it is accurate away from the smeared charges.

$\psi_\text{rest}$ is interpolated between corners with bicubic (Catmull–Rom) interpolation, so its
gradient (the field the lines follow) is continuous.

### 6.2 Causal field lines

$\psi$ is defined only up to a constant, and the constant decides where the lines are drawn. For a
static field it is chosen so that $\psi_\text{rest}$ has zero mean, which respects the
configuration's symmetries. But if the dynamic solver made the same choice every frame, dragging a
charge would shift *every* line on the canvas at once, including those the news of the motion
hasn't reached yet.

The fix uses the physics. Away from currents, the flux between two points changes as

$$\frac{\partial}{\partial t}\int_{\vec a}^{\vec b} \vec E\cdot\hat n\,\mathrm{d}l
= \int_{\vec a}^{\vec b} (\nabla\times B_z\hat z)\cdot\hat n\,\mathrm{d}l
= \pm\big[B_z(\vec b) - B_z(\vec a)\big],$$

so, with $\Phi = 2\pi K$ the flux per unit charge,

$$\Phi\,\frac{\partial\psi}{\partial t} = B_z + \text{const}.$$

**Fixing the constant.** Each frame, the constant is chosen so that $\psi$ at every grid corner
away from the charges has changed by $\int B_z\,\mathrm{d}t/\Phi$ since the previous frame.
- **Ahead of the wavefront** $B_z = 0$, so $\psi$ and the field lines there don't move at all.
  Causality becomes visible: the far-field lines only respond when the radiation reaches them.
- **Robust averaging.** The levels are periodic with period $1/N$, so the best constant is a
  circular mean over all corners.
- **Settling.** Once nothing has moved for a light-crossing time, the constant relaxes smoothly
  back to the static choice, so a settled field looks exactly like the static one.

### 6.3 Closed loops of radiation

A time-dependent field can have closed field lines: the familiar loops that detach from an
oscillating dipole. These touch no charge and needn't cross the canvas edge, so neither kind of
seed in section 3.4 finds them.

**Seeding.** Every closed contour of a smooth function encloses an extremum of it. So MEME finds
the extrema of $\psi$ on the grid and seeds along a ray from each extremum to the edge: every loop
round that extremum crosses the ray.

**Ill-conditioned extrema.** A long thin loop has an extremum on a nearly flat ridge, where
Newton's method can't pin it down. Saddles on the ridge, and grid corners where $\psi$ is a
discrete extremum among their eight neighbours, are used as seeds as well.

**Where the radiation hasn't arrived.** There the static fields of a neutral pair can cancel
almost exactly, and $\psi$ is numerical noise. Such regions (where $\lvert\nabla\psi\rvert$ is
negligible) are treated as having no field: no lines start there, and lines end on entering them.

## 7. Forces, torques and Gaussian surfaces

- **Forces on test charges.** A test charge feels $\vec F = q\vec E$ but makes no field. A test
  current $I$ in a field $\vec B$ in the plane feels $\vec F = I\,\hat\ell\times\vec B$ per unit
  length, with $\hat\ell$ along the wire.
- **Forces between wires.** These are the reverse of those between the equivalent charges: parallel
  currents attract.
- **Torques.** A rigid dipole (two charges or currents joined by `dipoleWith`) is drawn with the
  torque about its midpoint, $\tau = \sum_i (\vec r_i - \vec r_\text{mid}) \times \vec F_i$.
- **Gaussian surfaces.** A drawn closed curve shows an arrow at each crossing with a field line,
  outwards or inwards. In 2D, outward minus inward crossings equal $N$ times the enclosed charge
  exactly (section 3.2); for magnetic field lines they always cancel.

## 8. Where to find things

| Topic | Files |
|---|---|
| Analytic 2D models, normalisation | `fieldmodels.js`, `solvers.js` |
| Field lines, equipotentials, magnetic lines | `fieldlines2d.js` (seeding, bookkeeping), `levelsets2d.js` (tracing, stagnation points) |
| Grid, charge deposition, Maxwell solver | `maxwell/src/lib.rs`, `maxwell/src/stencil.rs` |
| Free-space static field | `maxwell/src/freespace.rs` |
| Absorbing layer | `maxwell/src/absorber.rs` |
| Stream function of a grid field, causal choice | `maxwell/src/streamfn.rs`, `gridfield.js` |
| Drawing, arrows, forces, Gaussian surfaces | `draw.js`, `fieldline.js`, `gaussian_surface.js` |
