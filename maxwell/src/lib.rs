use wasm_bindgen::prelude::*;
use serde::{Serialize, Deserialize};
#[allow(unused)]
use serde_wasm_bindgen::{from_value, to_value};
use console_error_panic_hook;
use web_sys::console as console;
use ndarray::Array2;
use std::fmt::{self, Debug, Formatter};
use crate::geometry::Geometry;

mod stencil;
mod fourier;
mod geometry;
mod contours;
mod streamfn;

const FIELD_SCALING: f64 = 20000.0;
const SOFTEN: f64 = 5.0;

#[wasm_bindgen]
pub struct Pair {
    pub u: f64,
    pub v: f64,
}


#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

#[wasm_bindgen]
#[derive(Serialize, Deserialize, Clone, PartialEq)]
pub struct Charge {
    x: f64,
    y: f64,
    charge: f64
}

#[wasm_bindgen]
impl Charge {
    #[wasm_bindgen(constructor)]
    pub fn new(x: f64, y: f64, charge: f64) -> Charge {
        Charge { x, y, charge }
    }
}


impl Debug for Charge {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "Charge {{ x: {}, y: {}, charge: {} }}", self.x, self.y, self.charge)
    }
}


#[wasm_bindgen]
pub struct FieldConfiguration {
    charges: Vec<Charge>,
    charges_at_last_tick: Vec<Charge>, 
    geometry: Geometry,
    cic_grid: Option<Array2<f64>>,
    elec_x: Option<Array2<f64>>,
    elec_y: Option<Array2<f64>>,
    current_x: Option<Array2<f64>>,
    current_y: Option<Array2<f64>>,
    mag_z: Option<Array2<f64>>,
    stencils: stencil::Stencils,
    charge_normalization: f64,
    elec_uniform_x: f64,
    elec_uniform_y: f64,
    stream_solver: Option<streamfn::StreamFunctionSolver>,
    /// running time-average of E, used by the absorbing layer (see tick)
    elec_avg_x: Option<Array2<f64>>,
    elec_avg_y: Option<Array2<f64>>,
    ticks_since_divergence_cleaning: usize,
    divergence_dirty: bool,
    currents_nonzero: bool,
    /// <E> after the last tick, and the corrections of its changes still to be applied, as
    /// (simulation time of the change, uniform change of E to apply) (see tick)
    mean_field: (f64, f64),
    pending_mean_corrections: std::collections::VecDeque<(f64, f64, f64)>,
    simulation_time: f64,
    // For the causal choice of field lines (streamfn.rs): time integrals of B_z and of the mean
    // current since psi was last computed, the time since then and since a charge last moved, and psi
    // at that time.
    b_integral: Option<Array2<f64>>,
    /// (uniform changes of E applied since psi was last computed)
    mean_current_integral: (f64, f64),
    time_since_stream_function: f64,
    time_since_motion: f64,
    gauge_history: Option<streamfn::GaugeHistory>,
}

/// Damping rates in the absorbing layer ramp linearly from 0 at its inner edge to this (per unit time).
const PML_MAX_SIGMA: f64 = 0.25;
/// Time constant of the running average of E in the absorbing layer: only E - <E> is damped, so
/// that static fields are not (internal time units; 200 = 1 s).
const PML_AVERAGE_TIME: f64 = 200.0;
/// Remove any violation of Gauss's law (created by the absorbing layer) every this many ticks.
const DIVERGENCE_CLEANING_TICKS: usize = 20;

pub fn evaluate_grid(field: &Array2<f64>, x: isize, y: isize) -> f64 {
    // Evaluate the field at the specified grid cell, or return 0 if x>=nx, y>=ny, x<0 or y<0
    let nx = field.shape()[0] as isize;
    let ny = field.shape()[1] as isize;
    if x < 0 || x >= nx || y < 0 || y >= ny {
        return 0.0;
    }
    field[[x as usize, y as usize]]
}

pub fn evaluate_grid_interpolated_or_0(field_config: &FieldConfiguration, grid: &Option<Array2<f64>>, x: f64, y: f64) -> f64 {
    if let Some(grid) = grid {
        evaluate_grid_interpolated(field_config, grid, x, y)
    } else {
        0.0
    }
}

pub fn evaluate_grid_interpolated(field_config: &FieldConfiguration, grid: &Array2<f64>, x: f64, y: f64) -> f64 {
    let nx = grid.shape()[0];
    let ny = grid.shape()[1];

    let dx = field_config.geometry.delta_x();
    let dy = field_config.geometry.delta_y();

    assert!(nx == field_config.geometry.nx && ny == field_config.geometry.ny, "Grid size mismatch");

    if x < 0.0 || x > field_config.geometry.x_max || y < 0.0 || y > field_config.geometry.y_max {
        console::log_2(
            &"Point out of bounds:".into(),
            &format!("({}, {})", x, y).into()
        );
        return 0.0;
    }

    let (i0, j0) = field_config.geometry.position_to_cell_unclamped(x, y);
    let (x0, y0) = field_config.geometry.cell_to_centroid(i0 as usize, j0 as usize);

    // if x>x0, we interpolate between i and i+1; if x<x0, we interpolate between i-1 and i
    let i1 = if x > x0 { i0 + 1 } else { i0 - 1 };
    let j1 = if y > y0 { j0 + 1 } else { j0 - 1 };

    // work out the weights for the interpolation
    let wx1 = (x - x0).abs() / dx;
    let wy1 = (y - y0).abs() / dy;
    let wx0 = 1.0 - wx1;
    let wy0 = 1.0 - wy1;

    // interpolate the field
    let mut field = 0.0;
    field += evaluate_grid(grid, i0, j0) * wx0 * wy0;
    field += evaluate_grid(grid, i1, j0) * wx1 * wy0;
    field += evaluate_grid(grid, i0, j1) * wx0 * wy1;
    field += evaluate_grid(grid, i1, j1) * wx1 * wy1;

    field
}



#[wasm_bindgen]
impl FieldConfiguration {
    #[wasm_bindgen(constructor)]
    pub fn new(x_max: f64, y_max: f64, nx: usize, ny: usize) -> FieldConfiguration {
        let geometry = Geometry{x_max: x_max, y_max: y_max, nx: nx, ny: ny, nboundary: nx/8};
        let geometry_clone = geometry.clone();
        let cell_area = geometry.delta_x() * geometry.delta_y();
        let charge_normalization = 4000.0 / cell_area;
        FieldConfiguration { charges: vec![], charges_at_last_tick: vec![], geometry, 
            cic_grid: None, elec_x: None, elec_y: None, mag_z: None, current_x: None, current_y: None,
            stencils: stencil::Stencils::new(geometry_clone), charge_normalization, 
            elec_uniform_x: 0.0, elec_uniform_y: 0.0, stream_solver: None,
            elec_avg_x: None, elec_avg_y: None, ticks_since_divergence_cleaning: 0, divergence_dirty: false, currents_nonzero: false,
            mean_field: (0.0, 0.0), pending_mean_corrections: std::collections::VecDeque::new(), simulation_time: 0.0, b_integral: None, mean_current_integral: (0.0, 0.0),
            time_since_stream_function: 0.0, time_since_motion: f64::INFINITY, gauge_history: None }
    }

    pub fn set_charges(&mut self, charges: JsValue) {
        self.charges = match from_value(charges) {
            Ok(charges) => charges,
            Err(err) => {
                console::log_1(&format!("Error deserializing charges: {:?}", err).into());
                vec![]
            }
        };
    }

    pub fn set_uniform_field(&mut self, x: f64, y: f64) {
        self.elec_uniform_x = x;
        self.elec_uniform_y = y;
    }

    pub fn get_uniform_field(&self) -> Pair {
        Pair { u: self.elec_uniform_x, v: self.elec_uniform_y }
    }

    pub fn reset_fields(&mut self) {
        self.cic_grid = None;
        self.elec_x = None;
        self.elec_y = None;
        self.mag_z = None;
        self.current_x = None;
        self.current_y = None;
        self.elec_avg_x = None;
        self.elec_avg_y = None;
        self.gauge_history = None;
    }


    pub fn make_cic_grid(&mut self) {
        // Each charge is a Gaussian (sigma = 2 cells) centred on its exact position, so that the
        // field changes smoothly as a charge is dragged
        let mut grid = Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny));
        let (dx, dy) = (self.geometry.delta_x(), self.geometry.delta_y());
        for charge in &self.charges {
            if !self.in_physical_region(charge) {
                console::log_1(&format!("Charge out of bounds: {:?}", charge).into());
                continue;
            }
            let (sx, wx) = self.stencils.gaussian_weights(charge.x, dx);
            let (sy, wy) = self.stencils.gaussian_weights(charge.y, dy);
            self.stencils.add_separable(&mut grid, sx, &wx, sy, &wy, charge.charge * self.charge_normalization);
        }
        self.cic_grid = Some(grid);
    }

    fn in_physical_region(&self, charge: &Charge) -> bool {
        charge.x >= 0.0 && charge.x <= self.geometry.x_max && charge.y >= 0.0 && charge.y <= self.geometry.y_max
    }

    pub fn initialize_on_constraints(&mut self) {
        self.make_cic_grid();

        let (elec_x, elec_y) = self.stencils.grad_del_squared_inv(self.cic_grid.as_ref().unwrap());
        self.elec_avg_x = Some(elec_x.clone());
        self.elec_avg_y = Some(elec_y.clone());
        self.elec_y = Some(elec_y);
        self.elec_x = Some(elec_x);
        self.ticks_since_divergence_cleaning = 0;
        self.mag_z = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));

        self.current_x = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
        self.current_y = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
        self.currents_nonzero = false;
        self.pending_mean_corrections.clear();
        {
            let n = (self.geometry.nx * self.geometry.ny) as f64;
            self.mean_field = (self.elec_x.as_ref().unwrap().sum() / n, self.elec_y.as_ref().unwrap().sum() / n);
        }

        // a new static field: the field lines start from the canonical choice
        self.b_integral = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
        self.mean_current_integral = (0.0, 0.0);
        self.time_since_stream_function = 0.0;
        self.time_since_motion = f64::INFINITY;
        self.gauge_history = None;
      
        self.charges_at_last_tick = self.charges.clone();
        
    }


    pub fn ensure_initialized(&mut self) {
        // Basic case: initialize if we have no fields yet
        if self.elec_y.is_none() {
            self.initialize_on_constraints();
        }

        // During dynamic evolution, the total number of charges must remain constant and the strength of each charge must remain constant.
        // Otherwise Maxwell's constraint equations are violated; we therefore need to reinitialize the fields.

        let mut constraint_violation = false;

        if self.charges_at_last_tick.len() != self.charges.len() {
            constraint_violation = true;
        }
        
        for (charge_earlier, charge_now) in self.charges_at_last_tick.iter().zip(self.charges.iter()) {
            if charge_earlier.charge != charge_now.charge {
                constraint_violation = true;                
            }
        }

        if constraint_violation {
            self.initialize_on_constraints();
        }
    }

    pub fn evaluate_cic_grid_interpolated(&mut self, x: f64, y: f64) -> f64 {
        evaluate_grid_interpolated_or_0(self, &self.cic_grid, x, y)
    }

    pub fn make_currents(&mut self, delta_t: f64) {
        // The current density of the charges' motion since the last tick, constructed so that the
        // continuity equation holds exactly on the grid (Esirkepov's method for the separable
        // Gaussian deposit of make_cic_grid): then Gauss's law, div E = rho, is preserved.
        //   rho1 - rho0 = q [ (Wx1 - Wx0) (Wy0 + Wy1) / 2 + (Wx0 + Wx1) / 2 (Wy1 - Wy0) ]
        // and the two terms are the (backward-difference) divergences of jx and jy.
        let (dx, dy) = (self.geometry.delta_x(), self.geometry.delta_y());
        let mut jx = self.current_x.take().unwrap();
        let mut jy = self.current_y.take().unwrap();
        let mut moved = false;
        if self.currents_nonzero {
            jx.fill(0.0);
            jy.fill(0.0);
        }

        if self.charges.len() != self.charges_at_last_tick.len() {
            panic!("Number of charges changed since last tick");
        }

        for (charge_earlier, charge_now) in self.charges_at_last_tick.iter().zip(self.charges.iter()) {
            if charge_earlier.charge != charge_now.charge {
                panic!("Strength of a charge changed since last tick");
            }
            if (charge_earlier.x == charge_now.x && charge_earlier.y == charge_now.y)
                || !self.in_physical_region(charge_earlier) || !self.in_physical_region(charge_now) {
                continue;
            }
            moved = true;
            let q = charge_now.charge * self.charge_normalization;
            let (sx0, wx0) = self.stencils.gaussian_weights(charge_earlier.x, dx);
            let (sx1, wx1) = self.stencils.gaussian_weights(charge_now.x, dx);
            let (sy0, wy0) = self.stencils.gaussian_weights(charge_earlier.y, dy);
            let (sy1, wy1) = self.stencils.gaussian_weights(charge_now.y, dy);
            // weights on a common index range
            let common = |s0: isize, w0: &[f64], s1: isize, w1: &[f64]| -> (isize, Vec<f64>, Vec<f64>) {
                let start = s0.min(s1);
                let end = (s0 + w0.len() as isize).max(s1 + w1.len() as isize);
                let mut a = vec![0.0; (end - start) as usize];
                let mut b = vec![0.0; (end - start) as usize];
                for (k, v) in w0.iter().enumerate() { a[(s0 - start) as usize + k] = *v; }
                for (k, v) in w1.iter().enumerate() { b[(s1 - start) as usize + k] = *v; }
                (start, a, b)
            };
            let (sx, ax0, ax1) = common(sx0, &wx0, sx1, &wx1);
            let (sy, ay0, ay1) = common(sy0, &wy0, sy1, &wy1);
            // jx[i] = -(dx / dt) q sum_{i' <= i} (Wx1 - Wx0)[i'] * (Wy0 + Wy1) / 2
            let mut cx = Vec::with_capacity(ax0.len());
            let mut acc = 0.0;
            for k in 0..ax0.len() { acc += ax1[k] - ax0[k]; cx.push(-acc * dx / delta_t); }
            let my: Vec<f64> = ay0.iter().zip(ay1.iter()).map(|(a, b)| 0.5 * (a + b)).collect();
            self.stencils.add_separable(&mut jx, sx, &cx, sy, &my, q);
            let mut cy = Vec::with_capacity(ay0.len());
            let mut acc = 0.0;
            for k in 0..ay0.len() { acc += ay1[k] - ay0[k]; cy.push(-acc * dy / delta_t); }
            let mx: Vec<f64> = ax0.iter().zip(ax1.iter()).map(|(a, b)| 0.5 * (a + b)).collect();
            self.stencils.add_separable(&mut jy, sx, &mx, sy, &cy, q);
        }

        self.currents_nonzero = moved;
        self.current_x = Some(jx);
        self.current_y = Some(jy);
        

    }
    
    /// Evolve the fields by one timestep
    pub fn tick(&mut self, delta_t: f64) {
        self.ensure_initialized();

        // compute the currents from any motion in the charges. 
        self.make_currents(delta_t);

        // save the charges so we can compute the currents next time
        self.charges_at_last_tick = self.charges.clone();

        // Now we are ready to evolve the fields by one timestep. Note that the B and density 
        // fields are half a tick behind and half a grid cell to the left of the E and j fields.
        let elec_x = self.elec_x.as_mut().unwrap();
        let elec_y = self.elec_y.as_mut().unwrap();
        let mag_z = self.mag_z.as_mut().unwrap();
        let avg_x = self.elec_avg_x.as_mut().unwrap();
        let avg_y = self.elec_avg_y.as_mut().unwrap();
        let current_x = self.current_x.as_ref().unwrap();
        let current_y = self.current_y.as_ref().unwrap();

        if (delta_t / self.geometry.delta_x()).abs() > 0.1 {
            console::log_1(&(format!("Warning: timestep is large compared to grid cell size: {}", delta_t / self.geometry.delta_x())).into());
        } 


        // The updates below are written over flat slices (index i * ny + j, periodic wrap at the
        // edges) for speed; they are the same as the stencils' forward/backward differences.
        let (nx, ny) = (self.geometry.nx, self.geometry.ny);
        let (dx_inv, dy_inv) = (1.0 / self.geometry.delta_x(), 1.0 / self.geometry.delta_y());
        // damping profiles of the absorbing layer (sigma_x depends on i only, sigma_y on j only)
        let profile = |n: usize| -> Vec<f64> {
            let w = self.geometry.nboundary;
            (0..n).map(|k| if k < w || k >= n - w { (1.0 - k.min(n - k - 1) as f64 / w as f64) * PML_MAX_SIGMA } else { 0.0 }).collect()
        };
        let (sig_x, sig_y) = (profile(nx), profile(ny));
        let ex = elec_x.as_slice_mut().unwrap();
        let ey = elec_y.as_slice_mut().unwrap();
        let bz = mag_z.as_slice_mut().unwrap();
        let ax = avg_x.as_slice_mut().unwrap();
        let ay = avg_y.as_slice_mut().unwrap();
        let jx = current_x.as_slice().unwrap();
        let jy = current_y.as_slice().unwrap();
        let ib = self.b_integral.as_mut().unwrap().as_slice_mut().unwrap();

        // Evolve Bz first, to get it from half a tick behind to half a tick ahead of E:
        // dB_z/dt = dE_x/dy - dE_y/dx (forward differences), then damp B in the absorbing layer
        // (static B is zero here: charges only ever move transiently)
        for i in 0..nx {
            let ip = if i + 1 < nx { i + 1 } else { 0 };
            for j in 0..ny {
                let jp = if j + 1 < ny { j + 1 } else { 0 };
                let k = i * ny + j;
                bz[k] += ((ex[i * ny + jp] - ex[k]) * dy_inv - (ey[ip * ny + j] - ey[k]) * dx_inv) * delta_t;
                let sigma = sig_x[i] + sig_y[j];
                if sigma > 0.0 { bz[k] *= 1.0 - sigma * delta_t; }
                ib[k] += bz[k] * delta_t;   // this B advances E below
            }
        }

        // Then E a whole tick (backward differences of B), and in the absorbing layer damp the
        // deviation of E from its running time-average, so that waves are absorbed but static fields
        // are untouched (damping E itself drains the Coulomb field of the charges through the boundary)
        let w_avg = delta_t / PML_AVERAGE_TIME;
        let mut damped: f64 = 0.0;
        for i in 0..nx {
            let im = if i > 0 { i - 1 } else { nx - 1 };
            for j in 0..ny {
                let jm = if j > 0 { j - 1 } else { ny - 1 };
                let k = i * ny + j;
                ex[k] += ((bz[k] - bz[i * ny + jm]) * dy_inv - jx[k]) * delta_t;
                ey[k] += (-(bz[k] - bz[im * ny + j]) * dx_inv - jy[k]) * delta_t;
                let (sx, sy) = (sig_x[i], sig_y[j]);
                if sx > 0.0 || sy > 0.0 {
                    ax[k] += (ex[k] - ax[k]) * w_avg;
                    ay[k] += (ey[k] - ay[k]) * w_avg;
                    let (ddx, ddy) = ((ex[k] - ax[k]) * sy * delta_t, (ey[k] - ay[k]) * sx * delta_t);
                    ex[k] -= ddx;
                    ey[k] -= ddy;
                    damped = damped.max(ddx.abs()).max(ddy.abs());
                }
            }
        }

        // Only damping that changed E can have violated Gauss's law. (Damping rates are at most
        // PML_MAX_SIGMA * delta_t ~ 0.2 per tick, so the explicit update above is stable.)
        if damped > 0.0 { self.divergence_dirty = true; }
        // On the periodic grid the mean of E is a global mode: moving a charge changes it (d<E>/dt =
        // -<j>; the absorbing layer changes it too), which in free space would not happen and would
        // leave a uniform field behind. Each change of <E> is undone by a uniform change of E, but only
        // after light has had time to cross the canvas: by then the wave that caused it has reached
        // everywhere, so nothing changes ahead of a wavefront and the settled field is the static one.
        self.simulation_time += delta_t;
        let n = (nx * ny) as f64;
        let mean_now = (self.elec_x.as_ref().unwrap().sum() / n, self.elec_y.as_ref().unwrap().sum() / n);
        let change = (mean_now.0 - self.mean_field.0, mean_now.1 - self.mean_field.1);
        if change != (0.0, 0.0) {
            self.pending_mean_corrections.push_back((self.simulation_time, -change.0, -change.1));
        }
        let delay = self.geometry.x_max.hypot(self.geometry.y_max);
        let (mut kick_x, mut kick_y) = (0.0, 0.0);
        while let Some(&(t, kx, ky)) = self.pending_mean_corrections.front() {
            if self.simulation_time - t < delay { break; }
            kick_x += kx; kick_y += ky;
            self.pending_mean_corrections.pop_front();
        }
        if kick_x != 0.0 || kick_y != 0.0 {
            // uniform: no curl, so no B; added to the layer's running average too, so not damped
            self.elec_x.as_mut().unwrap().mapv_inplace(|v| v + kick_x);
            self.elec_y.as_mut().unwrap().mapv_inplace(|v| v + kick_y);
            self.elec_avg_x.as_mut().unwrap().mapv_inplace(|v| v + kick_x);
            self.elec_avg_y.as_mut().unwrap().mapv_inplace(|v| v + kick_y);
            self.mean_current_integral.0 += kick_x;
            self.mean_current_integral.1 += kick_y;
        }
        self.mean_field = (mean_now.0 + kick_x, mean_now.1 + kick_y);
        self.time_since_stream_function += delta_t;
        self.time_since_motion = if self.currents_nonzero { 0.0 } else { self.time_since_motion + delta_t };

        self.ticks_since_divergence_cleaning += 1;
        if self.ticks_since_divergence_cleaning >= DIVERGENCE_CLEANING_TICKS {
            if self.divergence_dirty { self.clean_divergence(); }
            self.ticks_since_divergence_cleaning = 0;
        }
    }

    /// Restore Gauss's law, div E = rho (with backward differences, as the field update preserves it),
    /// by subtracting the gradient of del^-2 (div E - rho). The update conserves div E - rho exactly;
    /// only the damping in the absorbing layer violates it.
    pub fn clean_divergence(&mut self) {
        self.ticks_since_divergence_cleaning = 0;
        self.divergence_dirty = false;
        self.make_cic_grid();
        let elec_x = self.elec_x.as_ref().unwrap();
        let elec_y = self.elec_y.as_ref().unwrap();
        let rho = self.cic_grid.as_ref().unwrap();
        let mut error = Array2::<f64>::zeros(rho.dim());
        for i in 0..self.geometry.nx {
            for j in 0..self.geometry.ny {
                error[[i, j]] = self.stencils.evaluate(elec_x, i, j, &stencil::StencilType::GradX, &stencil::DifferenceType::Backward)
                    + self.stencils.evaluate(elec_y, i, j, &stencil::StencilType::GradY, &stencil::DifferenceType::Backward)
                    - rho[[i, j]];
            }
        }
        let (cx, cy) = self.stencils.grad_del_squared_inv(&error);
        *self.elec_x.as_mut().unwrap() -= &cx;
        *self.elec_y.as_mut().unwrap() -= &cy;
    }
}

impl FieldConfiguration {
    pub fn evaluate_elec_interpolated(&mut self, x: f64, y: f64) -> (f64, f64) {
        self.ensure_initialized();
        let elec_x = evaluate_grid_interpolated_or_0(self, &self.elec_x, x, y);
        let elec_y = evaluate_grid_interpolated_or_0(self, &self.elec_y, x, y);
        (elec_x + self.elec_uniform_x, elec_y + self.elec_uniform_y)
    }

    pub fn distance_to_closest_charge(&self, x: f64, y: f64) -> f64 {
        if(self.charges.len() == 0) {
            return 1.0e9;
        } else {
            let closest_charge = self.closest_charge(x, y).unwrap();
            let dx = x - closest_charge.x;
            let dy = y - closest_charge.y;
            (dx*dx + dy*dy).sqrt()
        }
    }

    pub fn closest_charge(&self, x: f64, y: f64) -> Option<&Charge> {
        let mut closest_charge: Option<&Charge> = None;
        let mut closest_distance: f64 = 1.0e9;
        for charge in &self.charges {
            let dx = x - charge.x;
            let dy = y - charge.y;
            let distance = (dx*dx + dy*dy).sqrt();
            if distance < closest_distance {
                closest_charge = Some(charge);
                closest_distance = distance;
            }
        }
        closest_charge
    }
}

#[wasm_bindgen]
pub fn compute_potential_electrostatic_direct(field_configuration: &FieldConfiguration, x: f64, y: f64) -> f64 {
    let mut potential: f64 = (x-field_configuration.geometry.x_max/2.0) * field_configuration.elec_uniform_x + (y - field_configuration.geometry.y_max/2.0) * field_configuration.elec_uniform_y;
    for charge in &field_configuration.charges {
        let dx = x - charge.x;
        let dy = y - charge.y;
        let r = (dx * dx + dy * dy + SOFTEN).sqrt();
        potential -= FIELD_SCALING * charge.charge / r;
    }
    potential
}

fn compute_field_electrostatic_direct_one_charge(charge: &Charge, x: f64, y: f64) -> Pair {
    let dx = x - charge.x;
    let dy = y - charge.y;
    let r = (dx * dx + dy * dy + SOFTEN).sqrt();
    let k = FIELD_SCALING * charge.charge;
    Pair { u: k * dx / (r * r * r), v: k * dy / (r * r * r) }
}

#[wasm_bindgen]
pub fn compute_field_electrostatic_direct(field_configuration: &FieldConfiguration, x: f64, y: f64) -> Pair {
    let mut u: f64 = field_configuration.elec_uniform_x;
    let mut v: f64 = field_configuration.elec_uniform_y;
    

    for charge in &field_configuration.charges {
        let Pair { u: u_i, v: v_i } = compute_field_electrostatic_direct_one_charge(charge, x, y);
        u += u_i;
        v += v_i;
    }
    Pair { u, v }
}

pub fn compute_force_electrostatic(field_config: &FieldConfiguration, for_charge: usize) -> Pair {
    let charge = &field_config.charges[for_charge];
    let mut u: f64 = 0.0;
    let mut v: f64 = 0.0;
    for other_charge in &field_config.charges {
        if other_charge as *const Charge == charge as *const Charge {
            continue;
        }
        let Pair { u: u_i, v: v_i } = compute_field_electrostatic_direct_one_charge(other_charge, charge.x, charge.y);
        u += u_i;
        v += v_i;
    }
    u*=charge.charge;
    v*=charge.charge;
    Pair { u, v }
}

#[wasm_bindgen]
pub fn compute_forces_electrostatic(field_config: &FieldConfiguration) -> Vec<Pair> {
    let mut forces = vec![];
    for i in 0..field_config.charges.len() {
        forces.push(compute_force_electrostatic(field_config, i));
    }
    forces
}

#[wasm_bindgen]
pub fn compute_one_force_electrostatic(field_config: &FieldConfiguration, x: f64, y: f64, charge: f64) -> Pair {
    let Pair { u, v } = compute_field_electrostatic_direct(field_config, x, y);
    Pair { u: u*charge, v: v*charge }
}

#[wasm_bindgen]
pub fn compute_one_force_magnetostatic(field_config: &FieldConfiguration, x: f64, y: f64, charge: f64) -> Pair {
    let Pair { u, v } = compute_field_magnetostatic_direct(field_config, x, y);
    Pair { u: v*charge, v: -u*charge }
}

#[wasm_bindgen]
pub fn compute_forces_magnetostatic(field_config: &FieldConfiguration) -> Vec<Pair> {
    let mut forces = compute_forces_electrostatic(field_config);
    // negate the electric force components to get the magnetic force in the plane
    for force in &mut forces {
        force.u = -force.u;
        force.v = -force.v;
    }
    forces
}

#[wasm_bindgen]
pub fn compute_field_electrostatic_direct_to_buffer(field_configuration: &FieldConfiguration, x: f64, y: f64, buffer: &mut [f64]) {
    let Pair { u, v } = compute_field_electrostatic_direct(field_configuration, x, y);
    buffer[0] = u;
    buffer[1] = v;
}

#[wasm_bindgen]
pub fn compute_field_electrostatic_per_charge_direct_to_buffer(field_configuration: &FieldConfiguration, x: f64, y: f64, buffer: &mut [f64]) {
    let mut offset = 0;
    for charge in 0..field_configuration.charges.len() {
        let Pair { u, v } = compute_field_electrostatic_direct_one_charge(&field_configuration.charges[charge], x, y);
        buffer[offset] = u;
        buffer[offset+1] = v;
        offset += 2;
    }
}

#[wasm_bindgen]
pub fn compute_field_magnetostatic_per_charge_direct_to_buffer(field_configuration: &FieldConfiguration, x: f64, y: f64, buffer: &mut [f64]) {
    let mut offset = 0;
    for charge in 0..field_configuration.charges.len() {
        let Pair { u, v } = compute_field_magnetostatic_direct_one_charge(&field_configuration.charges[charge], x, y);
        buffer[offset] = u;
        buffer[offset+1] = v;
        offset += 2;
    }
}



#[wasm_bindgen]
pub fn compute_field_magnetostatic_direct(field_configuration: &FieldConfiguration, x: f64, y:f64) -> Pair {
    // Here we take advantage of the fact the Biot-Savart law in a plane looks like $\hat{z} \times \nabla \phi$
    // where $\phi$ is a fictional electrostatic potential generated by replacing currents with charges.
    let Pair { u, v } = compute_field_electrostatic_direct(field_configuration, x, y);
    Pair { u: v, v: -u }
}

#[wasm_bindgen]
pub fn compute_field_magnetostatic_direct_one_charge(charge: &Charge, x: f64, y: f64) -> Pair {
    let Pair {u, v} = compute_field_electrostatic_direct_one_charge(charge, x, y);
    Pair { u: v, v: -u }
}

#[wasm_bindgen]
pub fn compute_field_magnetostatic_direct_to_buffer(field_configuration: &FieldConfiguration, x: f64, y: f64, buffer: &mut [f64]) {
    let Pair { u, v } = compute_field_magnetostatic_direct(field_configuration, x, y);
    buffer[0] = u;
    buffer[1] = v;
}

#[wasm_bindgen]
pub fn compute_electric_field_dynamic(field_configuration: &mut FieldConfiguration, x: f64, y: f64) -> Pair {
    let (u, v) = field_configuration.evaluate_elec_interpolated(x, y);
    Pair { u, v }
}

#[wasm_bindgen]
pub fn compute_electric_field_dynamic_to_buffer(field_configuration: &mut FieldConfiguration, x: f64, y: f64, buffer: &mut [f64]) {
    let Pair { u, v } = compute_electric_field_dynamic(field_configuration, x, y);
    buffer[0] = u;
    buffer[1] = v;
}

/// Single-valued part of the stream function of the grid electric field (see streamfn.rs), on the
/// corners of the physical region: (n+1)^2 values, x-major, where n = grid_physical_cells.
/// xs, ys, qs are the true positions and charges of the sources whose point-charge stream function
/// is added analytically. The additive constant chooses the drawn field lines: canonical for a
/// static field; for the dynamic solver it evolves causally from one call to the next (lines ahead
/// of the wavefront don't move), relaxing to the canonical choice once the field has settled.
#[wasm_bindgen]
pub fn compute_grid_stream_function_rest(field_configuration: &mut FieldConfiguration, xs: &[f64], ys: &[f64], qs: &[f64]) -> Vec<f64> {
    field_configuration.ensure_initialized();
    let fc = field_configuration;
    let solver = fc.stream_solver.get_or_insert_with(|| streamfn::StreamFunctionSolver::new(&fc.geometry));
    let flux_per_unit_charge = fc.charge_normalization * fc.geometry.delta_x() * fc.geometry.delta_y();
    let mut rest = streamfn::stream_function_rest(solver, &fc.geometry, fc.elec_x.as_ref().unwrap(), fc.elec_y.as_ref().unwrap(),
                                                  flux_per_unit_charge, xs, ys, qs);
    let mut raw = streamfn::free_psi_at_corners(&fc.geometry, xs, ys, qs);
    for (r, s) in raw.iter_mut().zip(rest.iter()) { *r += s; }
    let positions: Vec<(f64, f64)> = xs.iter().zip(ys.iter()).map(|(&x, &y)| (x, y)).collect();

    // the canonical choice: psi_rest as solved, with mean zero over the canvas
    let canonical = 0.0;
    let causal = match &fc.gauge_history {
        // only if psi was computed recently (every frame while field lines are shown)
        Some(history) if fc.time_since_stream_function <= GAUGE_MAX_GAP => {
            let deposits: Vec<(f64, f64)> = fc.charges_at_last_tick.iter().map(|c| (c.x, c.y)).collect();
            streamfn::causal_offset(&fc.geometry, &raw, history, fc.b_integral.as_ref().unwrap(), fc.mean_current_integral,
                                    flux_per_unit_charge, &positions, &deposits)
        }
        _ => None,
    };
    let offset = match causal {
        None => canonical,
        Some(c) => {
            // once the field has settled (no charge has moved for the time light takes to cross the
            // canvas), relax to the canonical choice, so that the result matches the static picture
            let settle_time = fc.geometry.x_max.hypot(fc.geometry.y_max);
            if fc.time_since_motion > settle_time {
                let w = 1.0 - (-fc.time_since_stream_function / GAUGE_RELAXATION_TIME).exp();
                c + w * streamfn::offset_difference(canonical, c)
            } else {
                c
            }
        }
    };

    for v in rest.iter_mut() { *v += offset; }
    for v in raw.iter_mut() { *v += offset; }
    fc.gauge_history = Some(streamfn::GaugeHistory { psi: raw, positions });
    fc.b_integral.as_mut().unwrap().fill(0.0);
    fc.mean_current_integral = (0.0, 0.0);
    fc.time_since_stream_function = 0.0;
    rest
}

/// Causal field lines need psi at least this often (internal time units; a frame is ~3.3)
const GAUGE_MAX_GAP: f64 = 20.0;
/// Time constant of the relaxation to the canonical field lines once the field has settled (0.5 s)
const GAUGE_RELAXATION_TIME: f64 = 100.0;

/// Number of grid cells across the physical region (the same in x and y).
#[wasm_bindgen]
pub fn grid_physical_cells(field_configuration: &FieldConfiguration) -> usize {
    field_configuration.geometry.nx - 2 * field_configuration.geometry.nboundary
}

#[cfg(test)]
mod stream_function_tests {
    use super::*;

    /// max over segments of |D psi - r| (outside `exclude`), absolute (in units of charge) and relative to max |r|
    fn misfit(fc: &FieldConfiguration, exclude: &dyn Fn(f64, f64) -> bool) -> (f64, f64) {
        let xs: Vec<f64> = fc.charges.iter().map(|c| c.x).collect();
        let ys: Vec<f64> = fc.charges.iter().map(|c| c.y).collect();
        let qs: Vec<f64> = fc.charges.iter().map(|c| c.charge).collect();
        let g = &fc.geometry;
        let flux = fc.charge_normalization * g.delta_x() * g.delta_y();
        let mut solver = streamfn::StreamFunctionSolver::new(g);
        let ex = fc.elec_x.as_ref().unwrap();
        let ey = fc.elec_y.as_ref().unwrap();
        let psi = streamfn::stream_function_rest(&mut solver, g, ex, ey, flux, &xs, &ys, &qs);
        let (rx, ry) = streamfn::residuals(g, ex, ey, flux, &xs, &ys, &qs);
        let m = g.nx - 2 * g.nboundary + 1;
        let (dx, dy) = (g.delta_x(), g.delta_y());
        let (mut worst, mut scale): (f64, f64) = (0.0, 1e-300);
        for a in 0..m {
            for b in 0..m {
                let p = a * m + b;
                if a + 1 < m {
                    scale = scale.max(rx[p].abs());
                    if !exclude((a as f64 + 0.5) * dx, b as f64 * dy) { worst = worst.max((psi[p + m] - psi[p] - rx[p]).abs()); }
                }
                if b + 1 < m {
                    scale = scale.max(ry[p].abs());
                    if !exclude(a as f64 * dx, (b as f64 + 0.5) * dy) { worst = worst.max((psi[p + 1] - psi[p] - ry[p]).abs()); }
                }
            }
        }
        (worst, worst / scale)
    }

    #[test]
    fn pure_radiation_is_reproduced_exactly() {
        // a coincident +-2 pair, one of them shaken: when they coincide again rho = 0 everywhere,
        // so the field is divergence-free and the stream function must reproduce every face flux
        let mut fc = FieldConfiguration::new(668.0, 335.0, 111, 111);
        let (x0, y0) = (334.0, 167.5);
        fc.charges = vec![Charge::new(x0, y0, 2.0), Charge::new(x0, y0, -2.0)];
        let step: f64 = 0.1 * 668.0 / 111.0;
        for frame in 1..=60 {
            let t = frame as f64 / 60.0;
            let y = if frame == 60 { y0 } else { y0 + 25.0 * (2.0 * std::f64::consts::PI * t / 0.5).sin() };
            fc.charges[0] = Charge::new(x0, y, 2.0);
            let mut dt: f64 = 1000.0 / 60.0 / 5.0;
            while dt > 0.0 { let s = step.min(dt); fc.tick(s); dt -= s; }
        }
        let (_, e) = misfit(&fc, &|_, _| false);
        assert!(e < 1e-9, "relative misfit {}", e);
    }

    #[test]
    fn psi_rest_is_continuous_as_a_charge_crosses_cells() {
        // drag a charge across cell boundaries in x and y in 0.25 px steps: psi_rest (and so every
        // field line) must change by a small amount at every step, with no jump at cell boundaries
        let mut prev: Option<Vec<f64>> = None;
        let mut steps = vec![];
        for k in 0..48 {
            let mut fc = FieldConfiguration::new(668.0, 335.0, 111, 111);
            fc.charges = vec![Charge::new(466.0 + 0.25 * k as f64, 165.0 + 0.125 * k as f64, 3.0), Charge::new(200.4, 167.5, -3.0)];
            fc.initialize_on_constraints();
            let xs: Vec<f64> = fc.charges.iter().map(|c| c.x).collect();
            let ys: Vec<f64> = fc.charges.iter().map(|c| c.y).collect();
            let qs: Vec<f64> = fc.charges.iter().map(|c| c.charge).collect();
            let flux = fc.charge_normalization * fc.geometry.delta_x() * fc.geometry.delta_y();
            let mut solver = streamfn::StreamFunctionSolver::new(&fc.geometry);
            let psi = streamfn::stream_function_rest(&mut solver, &fc.geometry, fc.elec_x.as_ref().unwrap(), fc.elec_y.as_ref().unwrap(), flux, &xs, &ys, &qs);
            if let Some(p) = prev {
                steps.push(psi.iter().zip(p.iter()).fold(0.0f64, |a, (u, v)| a.max((u - v).abs())));
            }
            prev = Some(psi);
        }
        let (lo, hi) = steps.iter().fold((f64::INFINITY, 0.0f64), |(a, b), v| (a.min(*v), b.max(*v)));
        // in units of the line spacing (1/4): the steps are all alike, and small
        assert!(4.0 * hi < 0.05 && hi < 3.0 * lo, "steps of psi_rest between {} and {} line spacings", 4.0 * lo, 4.0 * hi);
    }

    #[test]
    fn static_field_is_reproduced_away_from_the_charges() {
        // charges at cell centres (no snapping offset): outside the smearing kernels the misfit is
        // only the difference between the (anisotropic) smeared charge and a point charge, a small
        // fraction of the spacing of the field lines (1/4 in units of charge)
        let mut fc = FieldConfiguration::new(668.0, 335.0, 111, 111);
        let (dx, dy) = (fc.geometry.delta_x(), fc.geometry.delta_y());
        let c1 = (30.5 * dx, 40.5 * dy);
        let c2 = (55.5 * dx, 44.5 * dy);
        fc.charges = vec![Charge::new(c1.0, c1.1, 3.0), Charge::new(c2.0, c2.1, -2.0)];
        fc.initialize_on_constraints();
        let near = |x: f64, y: f64| [c1, c2].iter().any(|c| ((x - c.0) / dx).abs() < 12.0 && ((y - c.1) / dy).abs() < 12.0);
        let (abs, _) = misfit(&fc, &near);
        assert!(4.0 * abs < 0.05, "misfit {} line spacings", 4.0 * abs);
    }
}

#[cfg(test)]
mod dynamic_solver_tests {
    use super::*;

    /// Advance as ui.js does at 60 frames per second (ticks of 0.1 cell, internal time = ms / 5),
    /// calling `each_frame(fc, t_ms)` before each frame.
    fn run(fc: &mut FieldConfiguration, ms: f64, mut each_frame: impl FnMut(&mut FieldConfiguration, f64)) {
        let step: f64 = 0.1 * fc.geometry.x_max / fc.geometry.nx as f64;
        let frames = (ms / (1000.0 / 60.0)).round() as usize;
        for f in 1..=frames {
            each_frame(fc, f as f64 * 1000.0 / 60.0);
            let mut dt: f64 = 1000.0 / 60.0 / 5.0;
            while dt > 0.0 { let s = step.min(dt); fc.tick(s); dt -= s; }
        }
    }

    fn physical(fc: &FieldConfiguration) -> impl Iterator<Item = (usize, usize)> {
        let nb = fc.geometry.nboundary;
        let n = fc.geometry.nx - nb;
        (nb..n).flat_map(move |i| (nb..n).map(move |j| (i, j)))
    }

    /// max |E - E_poisson| over the physical region, relative to max |E_poisson|
    fn deviation_from_poisson(fc: &FieldConfiguration) -> f64 {
        let mut p = FieldConfiguration::new(fc.geometry.x_max, fc.geometry.y_max, fc.geometry.nx, fc.geometry.ny);
        p.charges = fc.charges.clone();
        p.initialize_on_constraints();
        let (ax, ay) = (fc.elec_x.as_ref().unwrap(), fc.elec_y.as_ref().unwrap());
        let (bx, by) = (p.elec_x.as_ref().unwrap(), p.elec_y.as_ref().unwrap());
        let (mut d, mut s): (f64, f64) = (0.0, 0.0);
        for (i, j) in physical(fc) {
            d = d.max((ax[[i, j]] - bx[[i, j]]).abs()).max((ay[[i, j]] - by[[i, j]]).abs());
            s = s.max(bx[[i, j]].abs()).max(by[[i, j]].abs());
        }
        d / s
    }

    fn pair() -> FieldConfiguration {
        let mut fc = FieldConfiguration::new(668.0, 335.0, 111, 111);
        fc.charges = vec![Charge::new(467.6, 167.5, 3.0), Charge::new(200.4, 167.5, -3.0)];
        fc.ensure_initialized();
        fc
    }

    #[test]
    fn static_charges_keep_their_field() {
        let mut fc = pair();
        run(&mut fc, 5000.0, |_, _| {});
        let d = deviation_from_poisson(&fc);
        assert!(d < 1e-9, "the static field drifted by {}", d);
    }

    #[test]
    fn after_a_drag_the_field_settles_to_the_static_field() {
        let mut fc = pair();
        run(&mut fc, 8000.0, |fc, t| {
            if t <= 300.0 { fc.charges[0] = Charge::new(467.6 + 60.0 * t / 300.0, 167.5 + 20.0 * t / 300.0, 3.0); }
        });
        let d = deviation_from_poisson(&fc);
        println!("deviation from the static field 7.7 s after the drag: {:.4}", d);
        assert!(d < 0.02, "deviation from the static field {}", d);
    }

    #[test]
    fn gauss_law_holds_during_a_sub_cell_drag() {
        // div E = rho (backward differences) in the physical region at every frame of a drag by
        // fractional amounts in both directions (the currents satisfy continuity exactly)
        let mut fc = pair();
        let mut worst: f64 = 0.0;
        run(&mut fc, 600.0, |fc, t| {
            fc.charges[0] = Charge::new(467.6 - 0.37 * t / 16.0, 167.5 + 0.23 * t / 16.0, 3.0);
            if t > 20.0 {
                let mut p = FieldConfiguration::new(668.0, 335.0, 111, 111);
                p.charges = fc.charges_at_last_tick.clone();
                p.make_cic_grid();
                let rho = p.cic_grid.as_ref().unwrap();
                let (ex, ey) = (fc.elec_x.as_ref().unwrap(), fc.elec_y.as_ref().unwrap());
                let (dx, dy) = (fc.geometry.delta_x(), fc.geometry.delta_y());
                let scale = rho.iter().fold(0.0f64, |a, v| a.max(v.abs()));
                for (i, j) in physical(fc) {
                    let div = (ex[[i, j]] - ex[[i - 1, j]]) / dx + (ey[[i, j]] - ey[[i, j - 1]]) / dy;
                    worst = worst.max((div - rho[[i, j]]).abs() / scale);
                }
            }
        });
        assert!(worst < 1e-9, "Gauss's law violated by {}", worst);
    }

    #[test]
    fn radiation_leaves_through_the_absorbing_layer() {
        // a coincident +-2 pair, one shaken for one period: rho = 0 afterwards, so all field is radiation
        let mut fc = FieldConfiguration::new(668.0, 335.0, 111, 111);
        fc.charges = vec![Charge::new(334.0, 167.5, 2.0), Charge::new(334.0, 167.5, -2.0)];
        fc.ensure_initialized();
        let energy = |fc: &FieldConfiguration| -> f64 {
            let (ex, ey, bz) = (fc.elec_x.as_ref().unwrap(), fc.elec_y.as_ref().unwrap(), fc.mag_z.as_ref().unwrap());
            physical(fc).map(|(i, j)| ex[[i, j]].powi(2) + ey[[i, j]].powi(2) + bz[[i, j]].powi(2)).sum()
        };
        let (mut peak, mut at4, mut at8): (f64, f64, f64) = (0.0, 0.0, 0.0);
        run(&mut fc, 8000.0, |fc, t| {
            let y = if t < 500.0 { 167.5 + 25.0 * (2.0 * std::f64::consts::PI * t / 500.0).sin() } else { 167.5 };
            fc.charges[0] = Charge::new(334.0, y, 2.0);
            let e = energy(fc);
            peak = peak.max(e);
            if (t - 4000.0).abs() < 8.0 { at4 = e; }
            at8 = e;
        });
        println!("radiation energy / peak: {:.1e} at 4 s, {:.1e} at 8 s", at4 / peak, at8 / peak);
        assert!(at4 / peak < 5e-3 && at8 / peak < 1e-4, "{} {}", at4 / peak, at8 / peak);
    }
}

#[cfg(test)]
mod causal_field_line_tests {
    use super::*;
    const LPU: f64 = streamfn::LINES_PER_UNIT_CHARGE;
    const CENTRE: (f64, f64) = (300.0, 167.5);
    const ARM: f64 = 20.0;

    fn dipole(angle: f64) -> Vec<Charge> {
        let (c, s) = (angle.cos(), angle.sin());
        vec![Charge::new(CENTRE.0 + ARM * c, CENTRE.1 + ARM * s, 2.0), Charge::new(CENTRE.0 - ARM * c, CENTRE.1 - ARM * s, -2.0)]
    }

    /// The drawn psi (in units of the line spacing) at the corners, as JS evaluates it
    fn drawn_psi(fc: &mut FieldConfiguration) -> Vec<f64> {
        let xs: Vec<f64> = fc.charges.iter().map(|c| c.x).collect();
        let ys: Vec<f64> = fc.charges.iter().map(|c| c.y).collect();
        let qs: Vec<f64> = fc.charges.iter().map(|c| c.charge).collect();
        let rest = compute_grid_stream_function_rest(fc, &xs, &ys, &qs);
        let free = streamfn::free_psi_at_corners(&fc.geometry, &xs, &ys, &qs);
        rest.iter().zip(free.iter()).map(|(r, f)| LPU * (r + f)).collect()
    }

    fn frame(fc: &mut FieldConfiguration) {
        let step: f64 = 0.1 * fc.geometry.x_max / fc.geometry.nx as f64;
        let mut dt: f64 = 1000.0 / 60.0 / 5.0;
        while dt > 0.0 { let s = step.min(dt); fc.tick(s); dt -= s; }
    }

    /// Largest change (in line spacings, modulo whole spacings) of psi between a and b over the corners
    /// at least `r` px from the centre (and 2 corners in from the edge)
    fn far_change(fc: &FieldConfiguration, a: &[f64], b: &[f64], r: f64) -> (f64, usize) {
        let m = fc.geometry.nx - 2 * fc.geometry.nboundary + 1;
        let (dx, dy) = (fc.geometry.delta_x(), fc.geometry.delta_y());
        let (mut worst, mut n) = (0.0f64, 0);
        for i in 2..m - 2 { for j in 2..m - 2 {
            if (i as f64 * dx - CENTRE.0).hypot(j as f64 * dy - CENTRE.1) < r { continue; }
            let d = a[i * m + j] - b[i * m + j];
            worst = worst.max((d - d.round()).abs());
            n += 1;
        } }
        (worst, n)
    }

    #[test]
    fn field_lines_outside_the_light_cone_do_not_move() {
        // spin a +-2 dipole one full turn in 0.6 s; c = 1 px per internal time unit (= 5 ms)
        let mut fc = FieldConfiguration::new(668.0, 335.0, 111, 111);
        fc.charges = dipole(0.0);
        fc.ensure_initialized();
        let before = drawn_psi(&mut fc);
        let mut checked = 0;
        for f in 1..=36 {
            let t_ms = f as f64 * 1000.0 / 60.0;
            fc.charges = dipole(2.0 * std::f64::consts::PI * (t_ms / 600.0).min(1.0));
            frame(&mut fc);
            let now = drawn_psi(&mut fc);
            // light-cone radius, plus the dipole arm and the smeared charge (4 sigma = 8 cells in x)
            let cone = t_ms / 5.0 + ARM + 10.0 * fc.geometry.delta_x();
            let (worst, n) = far_change(&fc, &before, &now, cone);
            if n > 100 {
                checked += 1;
                assert!(worst < 0.01, "frame {}: psi outside the light cone (r > {:.0} px) moved by {:.4} line spacings", f, cone, worst);
            }
        }
        assert!(checked >= 20, "too few frames checked: {}", checked);
    }

    #[test]
    fn settled_field_lines_match_the_static_ones() {
        let mut fc = FieldConfiguration::new(668.0, 335.0, 111, 111);
        fc.charges = dipole(0.0);
        fc.ensure_initialized();
        drawn_psi(&mut fc);
        let mut psi = vec![];
        for f in 1..=600 {   // spin 3/4 turn in 0.45 s, then 9.5 s at rest (the explorer animates for 10 s)
            let t_ms = f as f64 * 1000.0 / 60.0;
            fc.charges = dipole(1.5 * std::f64::consts::PI * (t_ms / 450.0).min(1.0));
            frame(&mut fc);
            psi = drawn_psi(&mut fc);
        }
        let mut fresh = FieldConfiguration::new(668.0, 335.0, 111, 111);
        fresh.charges = fc.charges.clone();
        fresh.ensure_initialized();
        let reference = drawn_psi(&mut fresh);
        let (worst, _) = far_change(&fc, &psi, &reference, 80.0);
        assert!(worst < 0.03, "settled field lines differ from the static ones by {:.4} line spacings", worst);
    }

    #[test]
    fn repeated_calls_do_not_change_the_lines() {
        let mut fc = FieldConfiguration::new(668.0, 335.0, 111, 111);
        fc.charges = dipole(0.3);
        fc.ensure_initialized();
        drawn_psi(&mut fc);
        fc.charges = dipole(0.5);
        frame(&mut fc);
        let a = drawn_psi(&mut fc);
        let b = drawn_psi(&mut fc);
        let (worst, _) = far_change(&fc, &a, &b, 0.0);
        assert!(worst < 1e-9, "{}", worst);
    }
}
