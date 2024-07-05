use wasm_bindgen::prelude::*;
use serde::{Serialize, Deserialize};
use serde_wasm_bindgen::{from_value, to_value};
use console_error_panic_hook;
use web_sys::console as console;
use ndarray::{Array2};
use std::{f32::consts::E, fmt::{self, Debug, Formatter}};
use num::integer::gcd;

mod stencil;
mod fourier;
mod pml;

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
    console::log_1(&"Hello from Rust, panic hook has been initialised!".into());
}

#[wasm_bindgen]
#[derive(Serialize, Deserialize, Clone)]
pub struct Charge {
    x: f64,
    y: f64,
    charge: f64,
}

#[wasm_bindgen]
impl Charge {
    #[wasm_bindgen(constructor)]
    pub fn new(x: f64, y: f64, charge: f64) -> Charge {
        Charge { x, y, charge }
    }
}

impl Charge {
    pub fn get_location_on_grid(&self, geometry: &Geometry) -> (usize, usize) {
        let (i,j) = geometry.position_to_cell(self.x, self.y).unwrap();
        (i, j)
    }
}


impl Debug for Charge {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "Charge {{ x: {}, y: {}, charge: {} }}", self.x, self.y, self.charge)
    }
}


#[derive(Clone)]
pub struct Geometry {
    /// maximum x extent of the physical region. Physical region runs from 0 to x_max; boundary cells extend further.
    x_max: f64, 
    /// maximum y extent of the physical region. Physical region runs from 0 to y_max; boundary cells extend further.
    y_max: f64, 
    /// number of grid cells in x direction, including boundary cells (=2*nboundary)
    nx: usize, 
    /// number of grid cells in y direction, including boundary cells (=2*nboundary)
    ny: usize, 
    /// number of boundary cells on each side
    nboundary: usize, 
}

impl Geometry {
    fn delta_x(&self) -> f64 {
        self.x_max / (self.nx - 2 * self.nboundary) as f64
    }

    fn delta_y(&self) -> f64 {
        self.y_max / (self.ny - 2 * self.nboundary) as f64
    }

    fn position_to_cell(&self, x: f64, y: f64) -> Option<(usize, usize)> {
        let (i, j) = self.position_to_cell_unclamped(x, y);
        if i < 0 || i >= self.nx as isize || j < 0 || j >= self.ny as isize {
            None
        } else {
            Some((i as usize, j as usize))
        }
    }

    fn position_to_cell_unclamped(&self, x: f64, y: f64) -> (isize, isize) {
        let i = (x / self.x_max * (self.nx - 2*self.nboundary) as f64) as isize + self.nboundary as isize;
        let j = (y / self.y_max * (self.ny - 2*self.nboundary) as f64) as isize + self.nboundary as isize;
        (i, j)
    }

    fn cell_to_centroid(&self, i: usize, j: usize) -> (f64, f64) {
        let x = (i as f64 - self.nboundary as f64 + 0.5) * self.delta_x();
        let y = (j as f64 - self.nboundary as f64 + 0.5) * self.delta_y();
        (x, y)
    }

    fn cell_to_corners(&self, i: usize, j: usize) -> Vec<(f64, f64)> {
        let x0 = (i as f64 - self.nboundary as f64) * self.delta_x();
        let x1 = (i as f64 - self.nboundary as f64 + 1.0) * self.delta_x();
        let y0 = (j as f64 - self.nboundary as f64) * self.delta_y();
        let y1 = (j as f64 - self.nboundary as f64 + 1.0) * self.delta_y();
        
        vec![(x0, y0), (x1, y0), (x0, y1), (x1, y1)]
    }

    fn x_extent_including_boundary(&self) -> f64 {
        self.x_max * self.nx as f64 / (self.nx - 2*self.nboundary) as f64
    }

    fn y_extent_including_boundary(&self) -> f64 {
        self.y_max * self.ny as f64 / (self.ny - 2*self.nboundary) as f64
    }

    fn in_padding_region(&self, x: f64, y: f64) -> bool {
        let min_x = -(self.nboundary as f64) * self.delta_x();
        let max_x = (self.nx as f64 - self.nboundary as f64) * self.delta_x();
        let min_y = -(self.nboundary as f64) * self.delta_y();
        let max_y = (self.ny as f64 - self.nboundary as f64) * self.delta_y();
        x>=min_x && x<=max_x && y>=min_y && y<=max_y
    }

}

#[wasm_bindgen]
pub struct FieldConfiguration {
    charges: Vec<Charge>,
    charges_at_last_tick: Vec<Charge>, 
    geometry: Geometry,
    cic_grid: Option<Array2<f64>>,
    Ex_grid: Option<Array2<f64>>,
    Ey_grid: Option<Array2<f64>>,
    jx_grid: Option<Array2<f64>>,
    jy_grid: Option<Array2<f64>>,
    Bz_grid: Option<Array2<f64>>,
    Bz_integral_grid: Option<Array2<f64>>, /// Bz field integrated over time, for use in the PML boundary conditions
    stencils: stencil::Stencils,
    charge_normalization: f64,
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
    #[wasm_bindgen(constructor)]
    pub fn new(x_max: f64, y_max: f64, nx: usize, ny: usize) -> FieldConfiguration {
        let geometry = Geometry{x_max: x_max, y_max: y_max, nx: nx, ny: ny, nboundary: nx/8};
        let geometry_clone = geometry.clone();
        let cell_area = geometry.delta_x() * geometry.delta_y();
        let charge_normalization = 4000.0 / cell_area;
        FieldConfiguration { charges: vec![], charges_at_last_tick: vec![], geometry, 
            cic_grid: None, Ex_grid: None, Ey_grid: None, Bz_grid: None, Bz_integral_grid: None, jx_grid: None, jy_grid: None,
            stencils: stencil::Stencils::new(geometry_clone), charge_normalization }
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

    pub fn reset_fields(&mut self) {
        self.cic_grid = None;
        self.Ex_grid = None;
        self.Ey_grid = None;
        self.Bz_grid = None;
        self.Bz_integral_grid = None;
        self.jx_grid = None;
        self.jy_grid = None;
    }


    pub fn make_cic_grid(&mut self) {
        self.cic_grid = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
        
        if let Some(ref mut grid) = self.cic_grid {
            for charge in &self.charges {
                if(charge.x<0.0 || charge.x>self.geometry.x_max || charge.y<0.0 || charge.y>self.geometry.y_max) {
                    console::log_1(&format!("Charge out of bounds: {:?}", charge).into());
                    continue;
                }
                let (i,j) = charge.get_location_on_grid(&self.geometry);
                self.stencils.add_softened_point(grid, i, j, charge.charge * self.charge_normalization);
            }
        }
    }

    pub fn initialize_on_constraints(&mut self) {
        self.make_cic_grid();

        let mut Ey: Array2<f64> = self.cic_grid.as_ref().unwrap().clone();
        let mut Ex = Ey.clone();
        self.stencils.apply(&mut Ey, stencil::StencilType::GradYDelSquaredInv, stencil::DifferenceType::Forward);
        self.stencils.apply(&mut Ex, stencil::StencilType::GradXDelSquaredInv, stencil::DifferenceType::Forward);

        self.Ey_grid = Some(Ey);
        self.Ex_grid = Some(Ex);
        self.Bz_grid = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
        self.Bz_integral_grid = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));

        self.jx_grid = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
        self.jy_grid = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
      
        self.charges_at_last_tick = self.charges.clone();
        
    }


    pub fn ensure_initialized(&mut self) {
        // Basic case: initialize if we have no fields yet
        if self.Ey_grid.is_none() {
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
        // Compute the current density from the charge density and the charge density at the last tick.
        // 
        // This is calculated by finding an approximate straight line between the locations of the charges
        
        let jx = self.jx_grid.as_mut().unwrap();
        let jy = self.jy_grid.as_mut().unwrap();
        
        jx.fill(0.0);
        jy.fill(0.0);
    
        if self.charges.len() != self.charges_at_last_tick.len() {
            panic!("Number of charges changed since last tick");
        }

        for (charge_earlier, charge_now) in self.charges_at_last_tick.iter().zip(self.charges.iter()) {

            if charge_earlier.charge != charge_now.charge {
                panic!("Strength of a charge changed since last tick");
            }

            let (i_earlier, j_earlier) = charge_earlier.get_location_on_grid(&self.geometry);
            let (i_now, j_now) = charge_now.get_location_on_grid(&self.geometry);

            if i_earlier == i_now && j_earlier == j_now {
                continue;
            }

            let delta_i = i_now as i32 - i_earlier as i32; // endpoint has di == delta_i
            let delta_j = j_now as i32 - j_earlier as i32; // endpoint has dj == delta_j

            let i_min = i_earlier.min(i_now);
            let i_max = i_earlier.max(i_now);
            let j_min = j_earlier.min(j_now);
            let j_max = j_earlier.max(j_now);

            // The current should be the charge density times the velocity averaged over the timestep
            // Note that if we are stepping more than one cell, the velocity is increased by a factor that is
            // precisely offset by the decrease in the weighting from time-averaging

            let x_current_density = (delta_i.signum() as f64) * charge_now.charge * self.charge_normalization * self.geometry.delta_x() / delta_t;
            let y_current_density = (delta_j.signum() as f64) * charge_now.charge * self.charge_normalization * self.geometry.delta_y() / delta_t;
            
            // Due to the relative grid alignment of current compared to charge fields, the current always runs from i_min inclusive 
            // to i_max exclusive (regardless of the direction of motion). Similarly for j_min and j_max.

            // On the assumption that the path length in a single frame is small, we just go along the x then the y direction rather
            // than anything fancier than that.

            for i in i_min..i_max {
                self.stencils.add_softened_point(jx, i, j_earlier, x_current_density);
            }

            for j in j_min..j_max {
                self.stencils.add_softened_point(jy, i_now, j, y_current_density);
            }
            

        }
        

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
        let Ex = self.Ex_grid.as_mut().unwrap();
        let Ey = self.Ey_grid.as_mut().unwrap();
        let Bz = self.Bz_grid.as_mut().unwrap();
        let Bz_integral = self.Bz_integral_grid.as_mut().unwrap();
        let jx = self.jx_grid.as_ref().unwrap();
        let jy = self.jy_grid.as_ref().unwrap();


        // Evolve Bz field first, to get it from half a tick behind to half a tick ahead of the E field
        for (i, j, sigma_x, sigma_y) in pml::pml_iterator_from_geometry(&self.geometry) {
            
            let d_Ex_dy = self.stencils.evaluate(&Ex, i, j, &stencil::StencilType::GradY, &stencil::DifferenceType::Forward);
            let d_Ey_dx = self.stencils.evaluate(&Ey, i, j, &stencil::StencilType::GradX, &stencil::DifferenceType::Forward);

            // dB_z/dt = -dE_x/dy + dE_y/dx
            Bz[[i,j]] += (d_Ex_dy - d_Ey_dx) * delta_t;
            
            // Now add the PML damping terms. This is a non-physical term that damps the fields near the boundary to 
            // mimic vacuum BCs. This is explained in the following sources:
            //  * https://onlinelibrary.wiley.com/doi/book/10.1002/9781118646700, chapter 3 (not very clearly written, but
            //    the equations are there)
            //  * https://arxiv.org/abs/2108.05348 (this is nicely clear but not explicit for the EM case)

            Bz[[i,j]] += (-(sigma_x+sigma_y) * Bz[[i, j]] - sigma_x*sigma_y * Bz_integral[[i,j]]) * delta_t;

            // Finally update the integral of Bz over time (which is just used for the PML damping, not anything physical)
            // Note that in tests I found neglecting the integral terms in the PML didn't make a huge difference to the
            // level of reflections (just a qualitative observation, not a rigorous test). Without the integral terms,
            // it's pretty obviously just a local damping term.
            
            Bz_integral[[i,j]] += Bz[[i,j]] * delta_t;
            

            
        }

        // Now the B field is half a tick ahead, so update the E field a tick to get ahead again
        for (i, j, sigma_x, sigma_y) in pml::pml_iterator_from_geometry(&self.geometry) {
            let d_Bz_dy = self.stencils.evaluate(&Bz, i, j, &stencil::StencilType::GradY, &stencil::DifferenceType::Backward);
            let d_Bz_dx = self.stencils.evaluate(&Bz, i, j, &stencil::StencilType::GradX, &stencil::DifferenceType::Backward);
            Ex[[i,j]] += (d_Bz_dy - jx[[i,j]]) * delta_t;
            Ey[[i,j]] += (-d_Bz_dx - jy[[i,j]]) * delta_t;
            
            // PML damping, as above:
            let pml_x_term = -sigma_y * Ex[[i, j]] + sigma_x * self.stencils.evaluate(&Bz_integral, i, j, &stencil::StencilType::GradY, &stencil::DifferenceType::Backward);
            let pml_y_term = -sigma_x * Ey[[i, j]] - sigma_y * self.stencils.evaluate(&Bz_integral, i, j, &stencil::StencilType::GradX, &stencil::DifferenceType::Backward) ;
            Ex[[i,j]] += pml_x_term*delta_t;
            Ey[[i,j]] += pml_y_term*delta_t;
        }

    }
}

impl FieldConfiguration {
    pub fn evaluate_E_grid_interpolated(&mut self, x: f64, y: f64) -> (f64, f64) {
        self.ensure_initialized();
        let Ex = evaluate_grid_interpolated_or_0(self, &self.Ex_grid, x, y);
        let Ey = evaluate_grid_interpolated_or_0(self, &self.Ey_grid, x, y);
        (Ex, Ey)
    }
}

#[wasm_bindgen]
pub fn compute_potential_electrostatic_direct(field_configuration: &FieldConfiguration, x: f64, y: f64) -> f64 {
    let mut potential: f64 = 0.0;
    for charge in &field_configuration.charges {
        let dx = x - charge.x;
        let dy = y - charge.y;
        let r = (dx * dx + dy * dy + SOFTEN).sqrt();
        potential -= FIELD_SCALING * charge.charge / r;
    }
    potential
}


pub fn generate_potential_contours(field_configuration: &FieldConfiguration, x0: f64, y0: f64, level: f64, reverse_direction: bool) -> Vec<(f64, f64)> {
    /// Generate a contour of the potential field starting at the specified point
    
    let STEP_SIZE: f64 = if reverse_direction { -0.5 } else { 0.5 };

    const RETURN_EVERY: i32 = 10;
    const FINISH_TOLERANCE: f64 = 0.5;
    const MAXIMUM_TRAVERSALS: i32 = 2;
    const MINIMUM_NUM_STEPS: i32 = 20;

    let mut x = x0;
    let mut y = y0;
    let mut contour: Vec<(f64, f64)> = vec![(x, y)];
    let mut step = 0;
    
    let mut needs_reversal = false; 

    loop {
        // predictor step...
        let Pair {u, v} = compute_field_electrostatic_direct(field_configuration, x, y);
        let r = (u * u + v * v).sqrt();
        
        x -= STEP_SIZE * v / r;
        y += STEP_SIZE * u / r;

        // corrector step...
        let Pair {u:u1, v: v1} = compute_field_electrostatic_direct(field_configuration, x, y);
        let r1 = (u1 * u1 + v1 * v1).sqrt();

        x -= STEP_SIZE * 0.5 * (v1 - v) / r1;
        y += STEP_SIZE * 0.5 * (u1 - u) / r1;

        step += 1;

        if step % RETURN_EVERY == 0 {
            contour.push((x, y));
        }
        if !field_configuration.geometry.in_padding_region(x, y) {
            if !reverse_direction {
                // If the contour doesn't close in on itself, we need to supplement it with the contour
                // running in the opposite direction from the start point
                let contour_reversed = generate_potential_contours(field_configuration, x0, y0, level, true);

                // now add contours_reversed onto the beginning of contour, in reverse order
                contour_reversed.iter().for_each(|(x, y)| contour.insert(0, (*x, *y))); 
            }

            break;
        }

        if ((x-x0).powi(2) + (y-y0).powi(2)) < FINISH_TOLERANCE.powi(2) && step > MINIMUM_NUM_STEPS {
            break;
        }
    }
    contour.push((x, y));
    contour
}


pub fn find_crossing_point(field_configuration: &FieldConfiguration, level: f64, x0: f64, y0: f64) -> (f64, f64) {
    /// Find a nearby point to (x0, y0) at which the potential has the specified level
    let mut x = x0;
    let mut y = y0;
    let mut potential = compute_potential_electrostatic_direct(field_configuration, x, y);
    let mut step = 0;

    const MAXIMUM_NUM_STEPS: i32 = 200;
    const TOLERANCE: f64 = 1e-4;
    
    let orig_potential = potential;
    

    if (potential-level).abs() < TOLERANCE {
        return (x0, y0);
    }
    let potential_offset_sign = (potential-level).signum();

    

    while (potential-level).signum() == potential_offset_sign {
        let Pair {u, v} = compute_field_electrostatic_direct(field_configuration, x, y);
        let r = (u * u + v * v).sqrt();
        let (u_normed, v_normed) = (u/r.powi(2), v/r.powi(2));

        if r<1e-6 || step>MAXIMUM_NUM_STEPS {
            console::log_1(&format!("Searching for crossing point from ({}, {}) at level {}; starting {}", x0, y0, level, orig_potential).into());
            console::log_1(&format!("Aborted after {} steps @ {} {} : {}", step, x, y, potential).into());

            return (x0, y0);
        }
        // try to overshoot the crossing point, then we'll bisect to find it exactly
        x -= 1.5*(potential-level)*u_normed;
        y -= 1.5*(potential-level)*v_normed;
        potential = compute_potential_electrostatic_direct(field_configuration, x, y);
        step += 1;
       
    } 
    let (mut x0, mut y0) = (x0, y0);
    let (mut x1, mut y1) = (x, y);

    let potential0 = compute_potential_electrostatic_direct(field_configuration, x0, y0);
    let potential1 = compute_potential_electrostatic_direct(field_configuration, x1, y1);

    // console::log_1(&format!("Found bisection range ({} {}):{} ({} {}):{} after {} steps", x0, y0, potential0, x1, y1, potential1, step).into());
    

    // now do bisection search between (x0, y0) and (x1, y1) to find crossing point
    step = 0;
   
    while (potential-level).abs() > TOLERANCE {
        x = (x0 + x1) / 2.0;
        y = (y0 + y1) / 2.0;
        potential = compute_potential_electrostatic_direct(field_configuration, x, y);

        if (potential-level).signum() == potential_offset_sign {
            x0 = x;
            y0 = y;
        } else {
            x1 = x;
            y1 = y;
        }
        
        step += 1;
        if step > MAXIMUM_NUM_STEPS {
            let final_pot = compute_potential_electrostatic_direct(field_configuration, x, y);
            console::log_1(&format!("Aborted bisection ({} {}): {} point after {} steps", x, y, final_pot, step).into());
            break;
        }
    }
    
    

    (x,y)
}

#[wasm_bindgen]
pub fn generate_potential_contours_at_level(field_configuration: &FieldConfiguration, level: f64) -> JsValue {
    /// Generate a contour at a specified level of the potential field
        
    const MAX_CONTOURS: usize = 5;
    
    // First, scan over the grid to find all cells that cross the level
    let mut crossing_flags = Array2::<bool>::from_shape_fn((field_configuration.geometry.nx, field_configuration.geometry.ny), |(i, j)| {
        let change_in_cell = field_configuration.geometry.cell_to_corners(i, j).
            iter().
            any(|(x, y)| compute_potential_electrostatic_direct(field_configuration, *x, *y) > level) !=
            field_configuration.geometry.cell_to_corners(i, j).
            iter().
            all(|(x, y)| compute_potential_electrostatic_direct(field_configuration, *x, *y) > level);

        change_in_cell
    });

    let mut contours : Vec<Vec<(f64, f64)>> = vec![];

    loop {
        if contours.len() >= MAX_CONTOURS {
            break;
        }


        // Find the first cell that crosses the level
        let (i, j) = match (0..field_configuration.geometry.nx).flat_map(|i| (0..field_configuration.geometry.ny).map(move |j| (i, j))).find(|(i, j)| crossing_flags[[*i, *j]]) {
            Some((i, j)) => (i, j),
            None => break,
        };

        crossing_flags[[i, j]] = false;

        // Now follow the contour from this cell
        let (mut x, mut y) = field_configuration.geometry.cell_to_centroid(i, j);
        (x, y) = find_crossing_point(field_configuration, level, x, y);
        
        
        let contour = generate_potential_contours(field_configuration, x, y, level, false);

        // unflag cells that have been visited by this contour
        contour.iter().for_each(|(x, y)| {
            match field_configuration.geometry.position_to_cell(*x, *y) {
                Some((it, jt)) => crossing_flags[[it, jt]] = false,
                None => (),
            }
        });


        contours.push(contour);


    }

    // web_sys::console::log_1(&format!("Found {} contours", contours.len()).into());

    to_value(&contours).unwrap()


}

#[wasm_bindgen]
pub fn compute_field_electrostatic_direct(field_configuration: &FieldConfiguration, x: f64, y: f64) -> Pair {
    let mut u: f64 = 0.0;
    let mut v: f64 = 0.0;
    

    for charge in &field_configuration.charges {
        let dx = x - charge.x;
        let dy = y - charge.y;
        let r = (dx * dx + dy * dy + SOFTEN).sqrt();
        let k = FIELD_SCALING * charge.charge;
        u += k * dx / (r * r * r);
        v += k * dy / (r * r * r);
    }
    Pair { u, v }
    // to_value(&(u, v)).unwrap()
}

#[wasm_bindgen]
pub fn compute_field_magnetostatic_direct(field_configuration: &FieldConfiguration, x: f64, y:f64) -> Pair {
    let mut u: f64 = 0.0;
    let mut v: f64 = 0.0;
    for current in &field_configuration.charges { // charges interpreted as current
        let dx = x - current.x;
        let dy = y - current.y;
        let r = (dx * dx + dy * dy).sqrt();
        let k = FIELD_SCALING * current.charge;
        u += k * dy / (r * r * r);
        v -= k * dx / (r * r * r);
    }
    Pair { u, v }
}

#[wasm_bindgen]
pub fn compute_electric_field_dynamic(field_configuration: &mut FieldConfiguration, x: f64, y: f64) -> Pair {
    let (u, v) = field_configuration.evaluate_E_grid_interpolated(x, y);
    Pair { u, v }
}