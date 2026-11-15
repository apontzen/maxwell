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
mod absorber;
mod freespace;

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
    /// the absorbing boundary layer (absorber.rs)
    absorber: Option<absorber::Absorber>,
    currents_nonzero: bool,
    simulation_time: f64,
    // For the causal choice of field lines (streamfn.rs): the time integral of B_z since psi was
    // last computed, the time since then and since a charge last moved, and psi at that time.
    b_integral: Option<Array2<f64>>,
    time_since_stream_function: f64,
    time_since_motion: f64,
    gauge_history: Option<streamfn::GaugeHistory>,
    // field lines per unit charge (the line spacing is 1 / this; see set_lines_per_unit_charge)
    lines_per_unit_charge: f64,
}


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
    /// Deprecated (cells are square only if nx / ny matches the region's aspect ratio): see
    /// with_cell_size. The physical region has (nx - 2 nx/8) x (ny - 2 nx/8) cells, surrounded by the
    /// absorbing layer.
    #[wasm_bindgen(constructor)]
    pub fn new(x_max: f64, y_max: f64, nx: usize, ny: usize) -> FieldConfiguration {
        FieldConfiguration::with_geometry(absorber::padded_geometry(x_max, y_max, nx, ny))
    }

    /// A grid of (nearly) square cells of side close to `cell` covering the region x_max x y_max
    /// (the nearest whole numbers of cells in x and y), surrounded by the absorbing layer.
    pub fn with_cell_size(x_max: f64, y_max: f64, cell: f64) -> FieldConfiguration {
        let n = |extent: f64| ((extent / cell).round() as usize).max(4);
        FieldConfiguration::with_geometry(absorber::physical_geometry(x_max, y_max, n(x_max), n(y_max)))
    }
}

impl FieldConfiguration {
    pub(crate) fn with_geometry(geometry: Geometry) -> FieldConfiguration {
        let geometry_clone = geometry.clone();
        let cell_area = geometry.delta_x() * geometry.delta_y();
        let charge_normalization = 4000.0 / cell_area;
        FieldConfiguration { charges: vec![], charges_at_last_tick: vec![], geometry, 
            cic_grid: None, elec_x: None, elec_y: None, mag_z: None, current_x: None, current_y: None,
            stencils: stencil::Stencils::new(geometry_clone), charge_normalization, 
            elec_uniform_x: 0.0, elec_uniform_y: 0.0, stream_solver: None,
            absorber: None, currents_nonzero: false,
            simulation_time: 0.0, b_integral: None,
            time_since_stream_function: 0.0, time_since_motion: f64::INFINITY, gauge_history: None,
            lines_per_unit_charge: streamfn::DEFAULT_LINES_PER_UNIT_CHARGE }
    }
}

#[wasm_bindgen]
impl FieldConfiguration {

    pub fn set_charges(&mut self, charges: JsValue) {
        self.charges = match from_value(charges) {
            Ok(charges) => charges,
            Err(err) => {
                console::log_1(&format!("Error deserializing charges: {:?}", err).into());
                vec![]
            }
        };
    }

    /// The number of field lines per unit charge drawn (compute_grid_stream_function_rest chooses
    /// psi's additive constant modulo the line spacing). A change restarts the causal choice.
    pub fn set_lines_per_unit_charge(&mut self, n: f64) {
        if n != self.lines_per_unit_charge {
            self.lines_per_unit_charge = n;
            self.gauge_history = None;
        }
    }

    /// The time step to use for tick(): 0.3 of the Courant stability limit of the grid (larger steps
    /// make the absorbing layer reflect more)
    pub fn stable_timestep(&self) -> f64 {
        let (dx, dy) = (self.geometry.delta_x(), self.geometry.delta_y());
        0.3 / (1.0 / (dx * dx) + 1.0 / (dy * dy)).sqrt()
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
        self.absorber = None;
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

    /// The charges as deposited on the grid: (x, y, charge * normalisation)
    fn deposit(&self) -> Vec<(f64, f64, f64)> {
        self.charges.iter().map(|c| (c.x, c.y, c.charge * self.charge_normalization)).collect()
    }

    fn in_physical_region(&self, charge: &Charge) -> bool {
        charge.x >= 0.0 && charge.x <= self.geometry.x_max && charge.y >= 0.0 && charge.y <= self.geometry.y_max
    }

    pub fn initialize_on_constraints(&mut self) {
        self.make_cic_grid();

        // the static field in free space (freespace.rs)
        let deposit = self.deposit();
        let potential = freespace::potential(&self.geometry, self.cic_grid.as_ref().unwrap(), &absorber::deposited_totals(&self.geometry, &deposit));
        let (elec_x, elec_y) = self.stencils.grad_forward(&potential);
        self.absorber = Some(absorber::Absorber::new(&self.geometry, &self.stencils, &deposit, &potential, self.simulation_time));
        self.elec_y = Some(elec_y);
        self.elec_x = Some(elec_x);
        self.mag_z = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));

        self.current_x = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
        self.current_y = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
        self.currents_nonzero = false;

        // a new static field: the field lines start from the canonical choice
        self.b_integral = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
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

        // compute the currents from any motion in the charges
        self.make_currents(delta_t);

        // save the charges so we can compute the currents next time
        self.charges_at_last_tick = self.charges.clone();

        if delta_t > 2.0 * self.stable_timestep() {
            console::log_1(&(format!("Warning: timestep {} exceeds the Courant limit {}", delta_t, 2.0 * self.stable_timestep())).into());
        }

        // the absorbing layer's static target (absorber.rs)
        self.simulation_time += delta_t;
        let deposit = self.deposit();
        let ab = self.absorber.as_mut().unwrap();
        ab.advance(&self.geometry, &self.stencils, self.simulation_time, &deposit);

        // The updates below are written over flat slices (index i * ny + j, periodic wrap at the
        // edges) for speed. B_z is half a tick behind E; it is advanced with forward differences of E,
        // then E a whole tick with backward differences of B_z, which preserves div E = rho.
        let (nx, ny) = (self.geometry.nx, self.geometry.ny);
        let (dx_inv, dy_inv) = (1.0 / self.geometry.delta_x(), 1.0 / self.geometry.delta_y());
        // semi-implicit (Crank-Nicolson) damping over a tick: f = (1 - s dt/2) / (1 + s dt/2), g = 1 / (1 + s dt/2)
        let factors = |sig: &Vec<f64>| -> (Vec<f64>, Vec<f64>) {
            (sig.iter().map(|s| (1.0 - 0.5 * s * delta_t) / (1.0 + 0.5 * s * delta_t)).collect(),
             sig.iter().map(|s| 1.0 / (1.0 + 0.5 * s * delta_t)).collect())
        };
        let (fx, gx) = factors(&ab.sigma_x);
        let (fy, gy) = factors(&ab.sigma_y);
        let (fxh, gxh) = factors(&ab.sigma_x_half);
        let (fyh, gyh) = factors(&ab.sigma_y_half);
        // cells touched by the layer (in either half-position)
        let in_x: Vec<bool> = (0..nx).map(|i| ab.sigma_x[i] > 0.0 || ab.sigma_x_half[i] > 0.0).collect();
        let in_y: Vec<bool> = (0..ny).map(|j| ab.sigma_y[j] > 0.0 || ab.sigma_y_half[j] > 0.0).collect();
        let ex = self.elec_x.as_mut().unwrap().as_slice_mut().unwrap();
        let ey = self.elec_y.as_mut().unwrap().as_slice_mut().unwrap();
        let bz = self.mag_z.as_mut().unwrap().as_slice_mut().unwrap();
        let bzx = ab.bzx.as_slice_mut().unwrap();
        let bzy = ab.bzy.as_slice_mut().unwrap();
        let tx = ab.target_x.as_slice().unwrap();
        let ty = ab.target_y.as_slice().unwrap();
        let jx = self.current_x.as_ref().unwrap().as_slice().unwrap();
        let jy = self.current_y.as_ref().unwrap().as_slice().unwrap();
        let ib = self.b_integral.as_mut().unwrap().as_slice_mut().unwrap();

        // dB_z/dt = dE_x/dy - dE_y/dx; in the layer B_z = B_zx + B_zy with dB_zx/dt = -dE_y/dx and
        // dB_zy/dt = dE_x/dy, damped, and driven by the deviation of E from the static target
        for i in 0..nx {
            let ip = if i + 1 < nx { i + 1 } else { 0 };
            for j in 0..ny {
                let jp = if j + 1 < ny { j + 1 } else { 0 };
                let k = i * ny + j;
                let (kx, ky) = (ip * ny + j, i * ny + jp);
                if in_x[i] || in_y[j] {
                    let d_ey = (ey[kx] - ty[kx]) - (ey[k] - ty[k]);
                    let d_ex = (ex[ky] - tx[ky]) - (ex[k] - tx[k]);
                    bzx[k] = fx[i] * bzx[k] - gx[i] * d_ey * dx_inv * delta_t;
                    bzy[k] = fy[j] * bzy[k] + gy[j] * d_ex * dy_inv * delta_t;
                    bz[k] = bzx[k] + bzy[k];
                } else {
                    // outside the layer B_z isn't split (its halves would grow without bound in a static field)
                    bz[k] += ((ex[ky] - ex[k]) * dy_inv - (ey[kx] - ey[k]) * dx_inv) * delta_t;
                }
                ib[k] += bz[k] * delta_t;   // this B advances E below
            }
        }

        // E: dE_x/dt = dB_z/dy - j_x, dE_y/dt = -dB_z/dx - j_y; in the layer the deviation from the
        // static target is damped (which violates Gauss's law there while it happens, but the layer
        // relaxes to the target, which obeys it; a global correction of div E would be acausal and
        // radiates noise)
        for i in 0..nx {
            let im = if i > 0 { i - 1 } else { nx - 1 };
            for j in 0..ny {
                let jm = if j > 0 { j - 1 } else { ny - 1 };
                let k = i * ny + j;
                let cx = ((bz[k] - bz[i * ny + jm]) * dy_inv - jx[k]) * delta_t;
                let cy = (-(bz[k] - bz[im * ny + j]) * dx_inv - jy[k]) * delta_t;
                if in_x[i] || in_y[j] {
                    let new_x = tx[k] + fyh[j] * (ex[k] - tx[k]) + gyh[j] * cx;
                    let new_y = ty[k] + fxh[i] * (ey[k] - ty[k]) + gxh[i] * cy;
                    ex[k] = new_x;
                    ey[k] = new_y;
                } else {
                    ex[k] += cx;
                    ey[k] += cy;
                }
            }
        }

        self.time_since_stream_function += delta_t;
        self.time_since_motion = if self.currents_nonzero { 0.0 } else { self.time_since_motion + delta_t };

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
/// corners of the physical region: (nx+1)(ny+1) values, x-major ((a, b) at a (ny+1) + b), where
/// [nx, ny] = grid_physical_cells.
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
            streamfn::causal_offset(&fc.geometry, &raw, history, fc.b_integral.as_ref().unwrap(),
                                    flux_per_unit_charge, &positions, &deposits, fc.lines_per_unit_charge)
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
                c + w * streamfn::offset_difference(canonical, c, fc.lines_per_unit_charge)
            } else {
                c
            }
        }
    };

    for v in rest.iter_mut() { *v += offset; }
    for v in raw.iter_mut() { *v += offset; }
    fc.gauge_history = Some(streamfn::GaugeHistory { psi: raw, positions });
    fc.b_integral.as_mut().unwrap().fill(0.0);
    fc.time_since_stream_function = 0.0;
    rest
}

/// Causal field lines need psi at least this often (internal time units; a frame is ~3.3)
const GAUGE_MAX_GAP: f64 = 20.0;
/// Time constant of the relaxation to the canonical field lines once the field has settled (0.5 s)
const GAUGE_RELAXATION_TIME: f64 = 100.0;

/// Numbers of grid cells across the physical region in x and y.
#[wasm_bindgen]
pub fn grid_physical_cells(field_configuration: &FieldConfiguration) -> Vec<u32> {
    let g = &field_configuration.geometry;
    vec![(g.nx - 2 * g.nboundary) as u32, (g.ny - 2 * g.nboundary) as u32]
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
        let (mx, my) = streamfn::corner_dims(g);
        let (dx, dy) = (g.delta_x(), g.delta_y());
        let (mut worst, mut scale): (f64, f64) = (0.0, 1e-300);
        for a in 0..mx {
            for b in 0..my {
                let p = a * my + b;
                if a + 1 < mx {
                    scale = scale.max(rx[p].abs());
                    if !exclude((a as f64 + 0.5) * dx, b as f64 * dy) { worst = worst.max((psi[p + my] - psi[p] - rx[p]).abs()); }
                }
                if b + 1 < my {
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
        let mut fc = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
        let (x0, y0) = (334.0, 167.5);
        fc.charges = vec![Charge::new(x0, y0, 2.0), Charge::new(x0, y0, -2.0)];
        let step: f64 = fc.stable_timestep();
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
            let mut fc = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
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
        let mut fc = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
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
        let step: f64 = fc.stable_timestep();
        let frames = (ms / (1000.0 / 60.0)).round() as usize;
        for f in 1..=frames {
            each_frame(fc, f as f64 * 1000.0 / 60.0);
            let mut dt: f64 = 1000.0 / 60.0 / 5.0;
            while dt > 0.0 { let s = step.min(dt); fc.tick(s); dt -= s; }
        }
    }

    fn physical(fc: &FieldConfiguration) -> impl Iterator<Item = (usize, usize)> {
        let nb = fc.geometry.nboundary;
        let (nx, ny) = (fc.geometry.nx - nb, fc.geometry.ny - nb);
        (nb..nx).flat_map(move |i| (nb..ny).map(move |j| (i, j)))
    }

    /// max |E - E_poisson| over the physical region, relative to max |E_poisson|
    fn deviation_from_poisson(fc: &FieldConfiguration) -> f64 {
        let mut p = FieldConfiguration::with_geometry(fc.geometry.clone());
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
        let mut fc = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
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
                let mut p = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
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
        let mut fc = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
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
        assert!(at4 / peak < 2e-6 && at8 / peak < 1e-8, "{} {}", at4 / peak, at8 / peak);
    }
}

#[cfg(test)]
mod causal_field_line_tests {
    use super::*;
    const LPU: f64 = streamfn::DEFAULT_LINES_PER_UNIT_CHARGE;
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
        let step: f64 = fc.stable_timestep();
        let mut dt: f64 = 1000.0 / 60.0 / 5.0;
        while dt > 0.0 { let s = step.min(dt); fc.tick(s); dt -= s; }
    }

    /// Largest change (in line spacings, modulo whole spacings) of psi between a and b over the corners
    /// at least `r` px from the centre (and 2 corners in from the edge)
    fn far_change(fc: &FieldConfiguration, a: &[f64], b: &[f64], r: f64) -> (f64, usize) {
        let (mx, my) = streamfn::corner_dims(&fc.geometry);
        let (dx, dy) = (fc.geometry.delta_x(), fc.geometry.delta_y());
        let (mut worst, mut n) = (0.0f64, 0);
        for i in 2..mx - 2 { for j in 2..my - 2 {
            if (i as f64 * dx - CENTRE.0).hypot(j as f64 * dy - CENTRE.1) < r { continue; }
            let d = a[i * my + j] - b[i * my + j];
            worst = worst.max((d - d.round()).abs());
            n += 1;
        } }
        (worst, n)
    }

    #[test]
    fn field_lines_outside_the_light_cone_do_not_move() {
        // spin a +-2 dipole one full turn in 0.6 s; c = 1 px per internal time unit (= 5 ms)
        let mut fc = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
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
        let mut fc = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
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
        let mut fresh = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
        fresh.charges = fc.charges.clone();
        fresh.ensure_initialized();
        let reference = drawn_psi(&mut fresh);
        let (worst, _) = far_change(&fc, &psi, &reference, 80.0);
        assert!(worst < 0.03, "settled field lines differ from the static ones by {:.4} line spacings", worst);
    }

    #[test]
    fn repeated_calls_do_not_change_the_lines() {
        let mut fc = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
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
/// The absorbing layer against vacuum: the same simulation on a grid 5x larger (whose walls are too
/// far away for reflections to come back in time), compared over the canvas.
#[cfg(test)]
mod vacuum_comparison {
    use super::*;
    fn run_case(fc: &mut FieldConfiguration, off: (f64, f64), co: (usize, usize), pos: &dyn Fn(f64) -> Vec<(f64, f64, f64)>, np: (usize, usize)) -> Vec<Vec<f64>> {
        fc.charges = pos(0.0).iter().map(|&(x, y, q)| Charge::new(x + off.0, y + off.1, q)).collect();
        fc.ensure_initialized();
        let nb = fc.geometry.nboundary;
        let mut snaps = vec![];
        for f in 1..=360 {
            let t = f as f64 * 1000.0 / 60.0;
            fc.charges = pos(t).iter().map(|&(x, y, q)| Charge::new(x + off.0, y + off.1, q)).collect();
            let mut dt: f64 = 1000.0 / 60.0 / 5.0;
            while dt > 0.0 { let s = fc.stable_timestep().min(dt); fc.tick(s); dt -= s; }
            if f % 60 == 0 {
                let (ex, ey, bz) = (fc.elec_x.as_ref().unwrap(), fc.elec_y.as_ref().unwrap(), fc.mag_z.as_ref().unwrap());
                let mut v = vec![];
                for i in 0..np.0 { for j in 0..np.1 { let (a, b) = (i + nb + co.0, j + nb + co.1); v.push(ex[[a, b]]); v.push(ey[[a, b]]); v.push(bz[[a, b]]); } }
                snaps.push(v);
            }
        }
        // minus the static field of the final configuration
        let mut p = FieldConfiguration::with_geometry(fc.geometry.clone());
        p.charges = fc.charges.clone(); p.initialize_on_constraints();
        let (ex, ey) = (p.elec_x.as_ref().unwrap(), p.elec_y.as_ref().unwrap());
        let mut st = vec![];
        for i in 0..np.0 { for j in 0..np.1 { let (a, b) = (i + nb + co.0, j + nb + co.1); st.push(ex[[a, b]]); st.push(ey[[a, b]]); st.push(0.0); } }
        snaps.iter().map(|v| v.iter().zip(st.iter()).map(|(a, b)| a - b).collect()).collect()
    }
    #[test]
    fn the_canvas_behaves_like_vacuum() {
        let np = (111usize, 56usize); let (dx, dy) = (668.0 / np.0 as f64, 335.0 / np.1 as f64);
        let rms = |v: &[f64]| (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt();
        let shake = |t: f64| -> Vec<(f64, f64, f64)> { let y = if t < 300.0 { 167.5 + 25.0 * (2.0 * std::f64::consts::PI * t / 300.0).sin() } else { 167.5 }; vec![(334.0, y, 2.0), (334.0, 167.5, -2.0)] };
        let drag = |t: f64| -> Vec<(f64, f64, f64)> { let u = (t / 300.0).min(1.0); vec![(467.6 + 60.0 * u, 167.5 + 20.0 * u, 3.0), (200.4, 167.5, -3.0)] };
        // largest allowed difference from vacuum (relative to the peak of the transient) at 1..6 s:
        // radiation is absorbed almost perfectly; after a drag the retarded target of the absorbing
        // layer (the static field of retarded positions, low-pass filtered) differs a little from
        // the true induction field until the field has settled. (With the periodic static field the
        // grid had before, its images made this 10x larger.)
        let limits = [("shake", [1e-5, 3e-4, 3e-4, 2e-4, 2e-4, 2e-4]), ("drag", [1e-5, 2e-3, 8e-3, 1.2e-2, 1e-2, 6e-3])];
        for ((name, limit), pos) in limits.iter().zip([&shake as &dyn Fn(f64) -> Vec<(f64, f64, f64)>, &drag]) {
            let mut big = FieldConfiguration::with_geometry(absorber::physical_geometry(5.0 * 668.0, 5.0 * 335.0, 5 * np.0, 5 * np.1));
            let r = run_case(&mut big, (2.0 * np.0 as f64 * dx, 2.0 * np.1 as f64 * dy), (2 * np.0, 2 * np.1), pos, np);
            let peak = r.iter().map(|v| rms(v)).fold(0.0f64, f64::max);
            let mut fc = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
            assert_eq!((fc.geometry.nx - 40, fc.geometry.ny - 40), np);
            let s = run_case(&mut fc, (0.0, 0.0), (0, 0), pos, np);
            for (k, (a, b)) in s.iter().zip(r.iter()).enumerate() {
                let e = rms(&a.iter().zip(b.iter()).map(|(x, y)| x - y).collect::<Vec<_>>()) / peak;
                println!("{} at {} s: difference from vacuum / peak {:.1e}", name, k + 1, e);
                assert!(e < limit[k], "{} at {} s: difference from vacuum {:.2e} (limit {:.0e})", name, k + 1, e, limit[k]);
            }
        }
    }
}

/// The grid's static field against the exact field of point charges in vacuum (2D: E = Q r / 2 pi r^2),
/// away from the charges (their Gaussian smearing and the grid's resolution matter within a few cells).
#[cfg(test)]
mod static_vacuum_tests {
    use super::*;

    /// Angles (degrees) between the grid's and the exact field, at the E_x and E_y nodes of the
    /// physical region further than `r_min` px from every charge
    fn angle_errors(charges: &[(f64, f64, f64)], r_min: f64) -> Vec<f64> {
        let mut fc = FieldConfiguration::with_cell_size(668.0, 335.0, 6.0);
        fc.charges = charges.iter().map(|&(x, y, q)| Charge::new(x, y, q)).collect();
        fc.reset_fields();
        fc.ensure_initialized();
        let g = fc.geometry.clone();
        let (dx, dy, nb) = (g.delta_x(), g.delta_y(), g.nboundary as f64);
        let (ex, ey) = (fc.elec_x.as_ref().unwrap(), fc.elec_y.as_ref().unwrap());
        let exact = |x: f64, y: f64| -> Option<(f64, f64)> {
            let (mut u, mut v) = (0.0, 0.0);
            for &(cx, cy, q) in charges {
                let (rx, ry) = (x - cx, y - cy);
                let r2 = rx * rx + ry * ry;
                if r2 < r_min * r_min { return None; }
                u += q * rx / r2; v += q * ry / r2;
            }
            Some((u, v))
        };
        let mut out = vec![];
        for i in g.nboundary..g.nx - g.nboundary - 1 {
            for j in g.nboundary..g.ny - g.nboundary - 1 {
                // the field at the corner (i + 1, j + 1): the mean of the neighbouring staggered components
                let (x, y) = ((i as f64 + 1.0 - nb) * dx, (j as f64 + 1.0 - nb) * dy);
                let Some((u, v)) = exact(x, y) else { continue };
                let gx = 0.5 * (ex[[i, j]] + ex[[i, j + 1]]);
                let gy = 0.5 * (ey[[i, j]] + ey[[i + 1, j]]);
                out.push((gx * v - gy * u).atan2(gx * u + gy * v).abs().to_degrees());
            }
        }
        out.sort_by(|a, b| a.partial_cmp(b).unwrap());
        out
    }

    #[test]
    fn the_static_field_is_that_of_charges_in_vacuum() {
        for (name, charges) in [
            ("horizontal dipole", vec![(0.3 * 668.0, 167.5, 2.0), (0.7 * 668.0, 167.5, -2.0)]),
            ("vertical dipole", vec![(334.0, 0.35 * 335.0, 2.0), (334.0, 0.65 * 335.0, -2.0)]),
            ("single charge", vec![(250.0, 120.0, 1.0)]),
        ] {
            let a = angle_errors(&charges, 40.0);
            let (median, p99, worst) = (a[a.len() / 2], a[a.len() * 99 / 100], a[a.len() - 1]);
            println!("{}: angle to the vacuum field beyond 40 px: median {:.2} deg, 99% {:.2}, max {:.2} ({} points)", name, median, p99, worst, a.len());
            assert!(median < 0.3 && p99 < 1.5, "{}: median {:.2} deg, 99% {:.2}", name, median, p99);
        }
    }
}
