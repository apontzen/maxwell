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
mod pml;
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
            elec_avg_x: None, elec_avg_y: None, ticks_since_divergence_cleaning: 0 }
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
        jx.fill(0.0);
        jy.fill(0.0);

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

        // Remove the mean current (a uniform neutralising background current, as is usual for periodic
        // domains). On the periodic grid the mean of E obeys d<E>/dt = -<j>, so otherwise moving a charge
        // would leave a uniform field behind; the divergence of j, and so Gauss's law, is unchanged.
        let n = jx.len() as f64;
        let (mean_x, mean_y) = (jx.sum() / n, jy.sum() / n);
        jx.mapv_inplace(|v| v - mean_x);
        jy.mapv_inplace(|v| v - mean_y);
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


        // Evolve Bz field first, to get it from half a tick behind to half a tick ahead of the E field
        for (i, j, sigma_x, sigma_y) in pml::pml_iterator_from_geometry(&self.geometry) {
            
            let d_elec_x_dy = self.stencils.evaluate(&elec_x, i, j, &stencil::StencilType::GradY, &stencil::DifferenceType::Forward);
            let d_elec_y_dx = self.stencils.evaluate(&elec_y, i, j, &stencil::StencilType::GradX, &stencil::DifferenceType::Forward);

            // dB_z/dt = -dE_x/dy + dE_y/dx
            mag_z[[i,j]] += (d_elec_x_dy - d_elec_y_dx) * delta_t;
            
            // Absorbing layer: damp B (static B is zero here: charges only ever move transiently)
            if sigma_x > 0.0 || sigma_y > 0.0 {
                mag_z[[i,j]] *= (-(sigma_x + sigma_y) * PML_MAX_SIGMA * delta_t).exp();
            }
        }

        // Now the B field is half a tick ahead, so update the E field a tick to get ahead again
        for (i, j, sigma_x, sigma_y) in pml::pml_iterator_from_geometry(&self.geometry) {
            let d_magz_dy = self.stencils.evaluate(&mag_z, i, j, &stencil::StencilType::GradY, &stencil::DifferenceType::Backward);
            let d_magz_dx = self.stencils.evaluate(&mag_z, i, j, &stencil::StencilType::GradX, &stencil::DifferenceType::Backward);
            elec_x[[i,j]] += (d_magz_dy - current_x[[i,j]]) * delta_t;
            elec_y[[i,j]] += (-d_magz_dx - current_y[[i,j]]) * delta_t;
            
            // Absorbing layer: damp the deviation of E from its running time-average, so that waves are
            // absorbed but static fields are untouched (damping E itself drains the Coulomb field of the
            // charges through the boundary)
            if sigma_x > 0.0 || sigma_y > 0.0 {
                let w = delta_t / PML_AVERAGE_TIME;
                avg_x[[i,j]] += (elec_x[[i,j]] - avg_x[[i,j]]) * w;
                avg_y[[i,j]] += (elec_y[[i,j]] - avg_y[[i,j]]) * w;
                elec_x[[i,j]] = avg_x[[i,j]] + (elec_x[[i,j]] - avg_x[[i,j]]) * (-sigma_y * PML_MAX_SIGMA * delta_t).exp();
                elec_y[[i,j]] = avg_y[[i,j]] + (elec_y[[i,j]] - avg_y[[i,j]]) * (-sigma_x * PML_MAX_SIGMA * delta_t).exp();
            }
        }

        self.ticks_since_divergence_cleaning += 1;
        if self.ticks_since_divergence_cleaning >= DIVERGENCE_CLEANING_TICKS {
            self.clean_divergence();
        }
    }

    /// Restore Gauss's law, div E = rho (with backward differences, as the field update preserves it),
    /// by subtracting the gradient of del^-2 (div E - rho). The update conserves div E - rho exactly;
    /// only the damping in the absorbing layer violates it.
    pub fn clean_divergence(&mut self) {
        self.ticks_since_divergence_cleaning = 0;
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
/// is added analytically.
#[wasm_bindgen]
pub fn compute_grid_stream_function_rest(field_configuration: &mut FieldConfiguration, xs: &[f64], ys: &[f64], qs: &[f64]) -> Vec<f64> {
    field_configuration.ensure_initialized();
    let fc = field_configuration;
    let solver = fc.stream_solver.get_or_insert_with(|| streamfn::StreamFunctionSolver::new(&fc.geometry));
    let flux_per_unit_charge = fc.charge_normalization * fc.geometry.delta_x() * fc.geometry.delta_y();
    streamfn::stream_function_rest(solver, &fc.geometry, fc.elec_x.as_ref().unwrap(), fc.elec_y.as_ref().unwrap(),
                                   flux_per_unit_charge, xs, ys, qs)
}

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
