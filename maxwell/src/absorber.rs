// The absorbing boundary of the dynamic solver: a perfectly matched layer (PML) round the physical
// region of the periodic grid. See https://github.com/apontzen/maxwell/issues/3.
//
// Berenger's split-field PML: in the layer B_z = B_zx + B_zy with
//     dB_zx/dt = -dE_y/dx - sigma_x B_zx,      dB_zy/dt = dE_x/dy - sigma_y B_zy,
//     dE_x/dt  =  dB_z/dy - j_x - sigma_y E_x,  dE_y/dt = -dB_z/dx - j_y - sigma_x E_y,
// which absorbs waves at any angle of incidence and frequency (in the continuum), with the
// conductivities graded as (depth / thickness)^3 and matched to the layer's thickness in x and y.
//
// A PML on its own also acts on static fields: in the layer a static E drives B_zx and B_zy to
// steady values -dE_y/dx / sigma_x and dE_x/dy / sigma_y that don't cancel where sigma_x != sigma_y,
// and the Coulomb field of the charges drains away. So here the layer acts only on the deviation of
// the fields from a static target E_t: E_x - E_t,x etc. (and B's static value, 0). The target is the
// static field of the charges where they were when the light now arriving at each cell left them
// (the retarded positions), so that the layer, like everything else, only learns that a charge
// has moved when the wave from it arrives. Once everything has settled the target is exactly the
// free-space static field of the charges (freespace.rs), which is then the steady state of the
// solver: inside the layer the field is pulled to the target, and a field that agrees with it there
// and has the same divergence (the charges) and curl (none) inside is the same field everywhere.
//
// The target is the (forward-difference) gradient of a retarded potential: the exact static
// potential of a reference configuration plus, for each charge that has moved recently,
// q * (G(cell - retarded position) - G(cell - reference position)), with G the static potential of a
// unit charge (precomputed, interpolated for sub-cell positions). Being a discrete gradient, the
// target has exactly zero discrete curl, so the total B_z = B_zx + B_zy still obeys Maxwell's
// equations exactly in the layer (only its halves see the target); a target mixing positions from
// different times would otherwise act as a source of B_z, radiating noise from the layer.
// When no charge has moved for longer than light takes to cross the grid, the reference is
// updated to the current positions (an exact static solve).

use std::collections::VecDeque;
use ndarray::Array2;
use crate::freespace;
use crate::geometry::Geometry;
use crate::stencil::Stencils;

/// Thickness of the layer, in cells.
pub const PML_CELLS: usize = 20;
/// Reflection coefficient the conductivity is designed for (at normal incidence, in the continuum)
const PML_REFLECTION: f64 = 1e-6;
/// Grading of the conductivity: sigma ~ depth^PML_GRADING
const PML_GRADING: i32 = 3;
/// Time constant of the low-pass filter on the target (internal time units; 200 = 1 s)
const TARGET_FILTER_TIME: f64 = 200.0;
/// The target is recomputed at least this often (internal time units; a frame is ~3.3)
const TARGET_INTERVAL: f64 = 2.0;

/// The geometry of a field configuration requested with an `nx` x `ny` grid: the physical region has
/// as many cells as it always has (nx - 2 nx/8, ny - 2 nx/8) and the layer PML_CELLS round it.
pub fn padded_geometry(x_max: f64, y_max: f64, nx: usize, ny: usize) -> Geometry {
    let old_boundary = nx / 8;
    physical_geometry(x_max, y_max, nx - 2 * old_boundary, ny - 2 * old_boundary)
}

/// The geometry with `npx` x `npy` cells in the physical region and the layer PML_CELLS round it.
pub fn physical_geometry(x_max: f64, y_max: f64, npx: usize, npy: usize) -> Geometry {
    Geometry { x_max, y_max, nx: npx + 2 * PML_CELLS, ny: npy + 2 * PML_CELLS, nboundary: PML_CELLS }
}

/// Conductivity profile along one axis with n cells, a layer of `w` cells each side of thickness
/// `thickness` px, evaluated at position k + offset (in cells; cell k spans [k, k + 1]), for k = 0..n.
/// Each field component is damped with the conductivity at its own (staggered) position.
fn profile(n: usize, w: usize, thickness: f64, offset: f64) -> Vec<f64> {
    let sigma_max = (PML_GRADING + 1) as f64 * (1.0 / PML_REFLECTION).ln() / (2.0 * thickness);
    let (wf, nf) = (w as f64, n as f64);
    (0..n).map(|k| {
        let x = k as f64 + offset;
        let depth = if x < wf { (wf - x) / wf } else if x > nf - wf { (x - (nf - wf)) / wf } else { 0.0 };
        sigma_max * depth.min(1.0).powi(PML_GRADING)
    }).collect()
}

/// The charges, (x, y, charge * charge_normalization), inside the physical region (the others are
/// not deposited: see FieldConfiguration::make_cic_grid), as (x, y, total charge of the deposit)
/// for freespace::potential.
pub fn deposited_totals(geometry: &Geometry, charges: &[(f64, f64, f64)]) -> Vec<(f64, f64, f64)> {
    let area = geometry.delta_x() * geometry.delta_y();
    charges.iter().filter(|&&(x, y, _)| x >= 0.0 && x <= geometry.x_max && y >= 0.0 && y <= geometry.y_max)
        .map(|&(x, y, q)| (x, y, q * area)).collect()
}

/// The free-space static potential (del^2 phi = rho, see freespace.rs) of point charges deposited as
/// in FieldConfiguration::make_cic_grid, whose forward-difference gradient is their static field.
/// Charges: (x, y, charge * charge_normalization).
pub fn static_potential(geometry: &Geometry, stencils: &Stencils, charges: &[(f64, f64, f64)]) -> Array2<f64> {
    let inside = deposited_totals(geometry, charges);
    let area = geometry.delta_x() * geometry.delta_y();
    let mut rho = Array2::<f64>::zeros((geometry.nx, geometry.ny));
    for &(x, y, q) in &inside {
        let (sx, wx) = stencils.gaussian_weights(x, geometry.delta_x());
        let (sy, wy) = stencils.gaussian_weights(y, geometry.delta_y());
        stencils.add_separable(&mut rho, sx, &wx, sy, &wy, q / area);
    }
    freespace::potential(geometry, &rho, &inside)
}

/// The free-space potential of a unit charge (in the units of `charges` of static_potential, i.e.
/// charge * charge_normalization = 1) at the centre of cell (nx, ny) of a grid twice the size
/// of `geometry` in each direction (so that shifting it by up to the grid's size never wraps), and
/// that cell.
fn unit_charge_table(geometry: &Geometry) -> (Array2<f64>, (f64, f64)) {
    let (dx, dy, w) = (geometry.delta_x(), geometry.delta_y(), geometry.nboundary);
    let (nx, ny) = (2 * geometry.nx, 2 * geometry.ny);
    let big = Geometry { x_max: (nx - 2 * w) as f64 * dx, y_max: (ny - 2 * w) as f64 * dy, nx, ny, nboundary: w };
    let stencils = Stencils::new(big.clone());
    let (gi, gj) = (geometry.nx, geometry.ny);
    let gx = (gi as f64 - w as f64 + 0.5) * dx;
    let gy = (gj as f64 - w as f64 + 0.5) * dy;
    (static_potential(&big, &stencils, &[(gx, gy, 1.0)]), (gi as f64, gj as f64))
}

/// A charge position prepared for lookups in the unit-charge table (see Absorber::prepare)
struct Prepared { di: usize, fa: f64, dj: usize, fb: f64 }

pub struct Absorber {
    /// conductivities at the positions of the field components (in cell units, cell (i, j) spanning
    /// [i, i+1] x [j, j+1]): E_x at (i+1, j+1/2), E_y at (i+1/2, j+1), B_z at (i+1, j+1). So E_x is
    /// damped with sigma_y_half[j], E_y with sigma_x_half[i], B_zx with sigma_x[i], B_zy with sigma_y[j]
    pub sigma_x: Vec<f64>,
    pub sigma_y: Vec<f64>,
    pub sigma_x_half: Vec<f64>,
    pub sigma_y_half: Vec<f64>,
    /// the split magnetic field in the layer (B_z = bzx + bzy there)
    pub bzx: Array2<f64>,
    pub bzy: Array2<f64>,
    /// the static target (meaningful in and next to the layer only), and its potential
    pub target_x: Array2<f64>,
    pub target_y: Array2<f64>,
    potential: Array2<f64>,
    /// the potential low-pass filtered in time (see advance), from which the target is taken
    potential_slow: Array2<f64>,
    pub(crate) filter_time: f64,
    /// cells where the potential is needed (within two cells of the layer), as (i, j), and those
    /// where the target is (within one)
    layer: Vec<(usize, usize)>,
    layer_target: Vec<(usize, usize)>,
    /// free-space potential of a unit charge at the centre of cell green_centre of a grid twice the
    /// size (unit_charge_table), made when first needed
    green: Option<Array2<f64>>,
    green_centre: (f64, f64),
    /// exact static potential of the reference configuration, and that configuration (x, y, q * norm)
    reference_potential: Array2<f64>,
    reference: Vec<(f64, f64, f64)>,
    /// for each charge, (time from which, x, y) of its positions over the last light-crossing time
    history: Vec<VecDeque<(f64, f64, f64)>>,
    last_motion: f64,
    last_update: f64,
    /// per charge and layer cell: the unit charge's field from its reference position, and from its
    /// current position (valid while the charge doesn't move: the time of its last move)
    green_reference: Vec<Vec<f64>>,
    green_current: Vec<(f64, Vec<f64>)>,
    /// the time light takes to cross the whole grid
    crossing_time: f64,
}

impl Absorber {
    /// charges: (x, y, charge * charge_normalization); static_potential: their exact static potential
    pub fn new(geometry: &Geometry, stencils: &Stencils, charges: &[(f64, f64, f64)], static_potential_: &Array2<f64>, time: f64) -> Absorber {
        let (nx, ny, w) = (geometry.nx, geometry.ny, geometry.nboundary);
        let (dx, dy) = (geometry.delta_x(), geometry.delta_y());
        let sigma_x = profile(nx, w, w as f64 * dx, 1.0);
        let sigma_y = profile(ny, w, w as f64 * dy, 1.0);
        let sigma_x_half = profile(nx, w, w as f64 * dx, 0.5);
        let sigma_y_half = profile(ny, w, w as f64 * dy, 0.5);
        // cells within `r` cells of the layer
        let near = |r: usize| -> Vec<(usize, usize)> {
            (0..nx).flat_map(|i| (0..ny).map(move |j| (i, j)))
                .filter(|&(i, j)| i < w + r || i + w + r >= nx || j < w + r || j + w + r >= ny)
                .collect()
        };
        let layer = near(2);
        let layer_target = near(1);
        let (target_x, target_y) = stencils.grad_forward(static_potential_);
        Absorber {
            sigma_x, sigma_y, sigma_x_half, sigma_y_half,
            bzx: Array2::zeros((nx, ny)), bzy: Array2::zeros((nx, ny)),
            target_x, target_y, potential: static_potential_.clone(), potential_slow: static_potential_.clone(), filter_time: TARGET_FILTER_TIME,
            layer, layer_target, green: None, green_centre: (0.0, 0.0),
            reference_potential: static_potential_.clone(), reference: charges.to_vec(),
            history: charges.iter().map(|&(x, y, _)| VecDeque::from(vec![(f64::NEG_INFINITY, x, y)])).collect(),
            last_motion: f64::NEG_INFINITY, last_update: time,
            green_reference: vec![vec![]; charges.len()],
            green_current: vec![(f64::NAN, vec![]); charges.len()],
            crossing_time: (nx as f64 * dx).hypot(ny as f64 * dy),
        }
    }

    /// For a charge at (x, y) in the physical region: the shift into the unit-charge table (integer
    /// and fractional parts), the same for every cell
    fn prepare(&self, geometry: &Geometry, x: f64, y: f64) -> Prepared {
        let w = geometry.nboundary as f64;
        // the charge's position in units of cells, relative to the centre of cell 0
        let u = x / geometry.delta_x() + w - 0.5;
        let v = y / geometry.delta_y() + w - 0.5;
        // (for a charge in the physical region and any cell of the grid, the shifted index and the
        // next one lie within the table)
        let (a, b) = ((self.green_centre.0 - u).max(0.0), (self.green_centre.1 - v).max(0.0));
        let (a0, b0) = (a.floor(), b.floor());
        Prepared { di: a0 as usize, fa: a - a0, dj: b0 as usize, fb: b - b0 }
    }

    /// The unit charge's potential at cell (i, j) for a charge at a prepared position: the table
    /// shifted by a (fractional) number of cells, interpolated bilinearly.
    #[inline]
    fn green_at(&self, i: usize, j: usize, p: &Prepared) -> f64 {
        let t = self.green.as_ref().expect("unit-charge table");
        let (nx, ny) = t.dim();
        let i0 = (i + p.di).min(nx - 2);
        let j0 = (j + p.dj).min(ny - 2);
        let (fa, fb) = (p.fa, p.fb);
        (1.0 - fa) * ((1.0 - fb) * t[[i0, j0]] + fb * t[[i0, j0 + 1]]) + fa * ((1.0 - fb) * t[[i0 + 1, j0]] + fb * t[[i0 + 1, j0 + 1]])
    }

    fn green(&self, geometry: &Geometry, i: usize, j: usize, x: f64, y: f64) -> f64 {
        let p = self.prepare(geometry, x, y);
        self.green_at(i, j, &p)
    }

    /// Low-pass filter (first order, time constant filter_time) of the retarded potential: the layer
    /// should only enforce the static part of the field, not the induction and radiation fields,
    /// which the retarded Coulomb potential doesn't describe
    fn filter(&mut self, elapsed: f64) {
        let w = if self.filter_time > 0.0 { 1.0 - (-elapsed / self.filter_time).exp() } else { 1.0 };
        for &(i, j) in &self.layer {
            let p = self.potential[[i, j]];
            self.potential_slow[[i, j]] += (p - self.potential_slow[[i, j]]) * w;
        }
    }

    /// The target from the potential: its forward-difference gradient, next to the layer
    fn target_from_potential(&mut self, geometry: &Geometry) {
        let (nx, ny) = (geometry.nx, geometry.ny);
        let (dx_inv, dy_inv) = (1.0 / geometry.delta_x(), 1.0 / geometry.delta_y());
        for &(i, j) in &self.layer_target {
            let (ip, jp) = ((i + 1) % nx, (j + 1) % ny);
            self.target_x[[i, j]] = (self.potential_slow[[ip, j]] - self.potential_slow[[i, j]]) * dx_inv;
            self.target_y[[i, j]] = (self.potential_slow[[i, jp]] - self.potential_slow[[i, j]]) * dy_inv;
        }
    }

    /// Advance to `time` with the charges (the grid's deposit) at `charges` = (x, y, q * norm).
    pub fn advance(&mut self, geometry: &Geometry, stencils: &Stencils, time: f64, charges: &[(f64, f64, f64)]) {
        let mut moved = false;
        for (k, &(x, y, _)) in charges.iter().enumerate() {
            let last = *self.history[k].back().unwrap();
            if (last.1, last.2) != (x, y) {
                self.history[k].push_back((time, x, y));
                moved = true;
            }
            // forget positions older than a crossing time (keeping the one in force then)
            while self.history[k].len() > 2 && self.history[k][1].0 < time - self.crossing_time { self.history[k].pop_front(); }
        }
        if moved { self.last_motion = time; }
        if !moved && time - self.last_update < TARGET_INTERVAL { return; }
        let elapsed = time - self.last_update;
        self.last_update = time;

        // once the waves of the last move have crossed the grid, the retarded potential is the static
        // potential of the charges where they are: use it exactly; and once the filter has caught up
        // with it, the target is exactly the static field
        let settled = time - self.last_motion > self.crossing_time;
        if settled && self.reference.iter().zip(charges.iter()).any(|(a, b)| (a.0, a.1) != (b.0, b.1)) {
            self.reference_potential = static_potential(geometry, stencils, charges);
            self.reference = charges.to_vec();
            for g in self.green_reference.iter_mut() { g.clear(); }
            for &(i, j) in &self.layer { self.potential[[i, j]] = self.reference_potential[[i, j]]; }
        }
        if settled {
            if time - self.last_motion > self.crossing_time + 8.0 * self.filter_time {
                if self.potential_slow != self.potential {
                    for &(i, j) in &self.layer { self.potential_slow[[i, j]] = self.potential[[i, j]]; }
                    self.target_from_potential(geometry);
                }
            } else {
                self.filter(elapsed);
                self.target_from_potential(geometry);
            }
            return;
        }

        if self.green.is_none() {
            let (table, centre) = unit_charge_table(geometry);
            self.green = Some(table);
            self.green_centre = centre;
        }
        // the unit charge's field from the reference positions, and from the current positions of
        // charges that have moved (cached until they move again)
        let n = self.layer.len();
        for k in 0..charges.len() {
            let (rx, ry, _) = self.reference[k];
            let moved_since_reference = self.history[k].iter().any(|e| (e.1, e.2) != (rx, ry));
            if moved_since_reference && self.green_reference[k].is_empty() {
                self.green_reference[k] = (0..n).map(|idx| { let (i, j) = self.layer[idx]; self.green(geometry, i, j, rx, ry) }).collect();
            }
            let last_move = self.history[k].back().unwrap().0;
            if moved_since_reference && !(self.green_current[k].0 == last_move) {
                let (x, y, _) = charges[k];
                let g = (0..n).map(|idx| { let (i, j) = self.layer[idx]; self.green(geometry, i, j, x, y) }).collect();
                self.green_current[k] = (last_move, g);
            }
        }

        // each charge's past positions, prepared
        let prepared: Vec<Vec<Prepared>> = (0..charges.len()).map(|k| {
            if self.green_reference[k].is_empty() { vec![] } else { self.history[k].iter().map(|e| self.prepare(geometry, e.1, e.2)).collect() }
        }).collect();
        let (dx, dy, w) = (geometry.delta_x(), geometry.delta_y(), geometry.nboundary as f64);
        for idx in 0..n {
            let (i, j) = self.layer[idx];
            let (cx, cy) = ((i as f64 - w + 0.5) * dx, (j as f64 - w + 0.5) * dy);
            let mut phi = self.reference_potential[[i, j]];
            for k in 0..charges.len() {
                if self.green_reference[k].is_empty() { continue; }   // hasn't moved since the reference
                let (x, y, q) = charges[k];
                let (rx, ry, _) = self.reference[k];
                // retarded time (c = 1 px per time unit), from the distance to where the charge is now
                let t_ret = time - (cx - x).hypot(cy - y);
                let last_move = self.history[k].back().unwrap().0;
                let a = if t_ret >= last_move {
                    self.green_current[k].1[idx]          // the wave of its last move has passed
                } else {
                    // between two recorded positions, interpolate linearly in time (a charge moved
                    // from one to the next over that interval), so that the target is continuous
                    let h = &self.history[k];
                    let e = h.partition_point(|e| e.0 <= t_ret).max(1) - 1;
                    if (h[e].1, h[e].2) == (rx, ry) && (e + 1 >= h.len() || (h[e + 1].1, h[e + 1].2) == (rx, ry)) { continue; }
                    let g0 = self.green_at(i, j, &prepared[k][e]);
                    if e + 1 < h.len() && h[e].0.is_finite() {
                        let f = ((t_ret - h[e].0) / (h[e + 1].0 - h[e].0)).clamp(0.0, 1.0);
                        g0 + f * (self.green_at(i, j, &prepared[k][e + 1]) - g0)
                    } else {
                        g0
                    }
                };
                phi += q * (a - self.green_reference[k][idx]);
            }
            self.potential[[i, j]] = phi;
        }
        self.filter(elapsed);
        self.target_from_potential(geometry);
    }
}
