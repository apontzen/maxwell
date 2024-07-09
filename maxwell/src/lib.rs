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
#[derive(Serialize, Deserialize, Clone, PartialEq)]
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
    mag_z_integral: Option<Array2<f64>>, /// Bz field integrated over time, for use in the PML boundary conditions
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
            cic_grid: None, elec_x: None, elec_y: None, mag_z: None, mag_z_integral: None, current_x: None, current_y: None,
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
        self.elec_x = None;
        self.elec_y = None;
        self.mag_z = None;
        self.mag_z_integral = None;
        self.current_x = None;
        self.current_y = None;
    }


    pub fn make_cic_grid(&mut self) {
        self.cic_grid = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
        
        if let Some(ref mut grid) = self.cic_grid {
            for charge in &self.charges {
                if charge.x<0.0 || charge.x>self.geometry.x_max || charge.y<0.0 || charge.y>self.geometry.y_max {
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

        let mut elec_y: Array2<f64> = self.cic_grid.as_ref().unwrap().clone();
        let mut elec_x = elec_y.clone();
        self.stencils.apply(&mut elec_y, stencil::StencilType::GradYDelSquaredInv, stencil::DifferenceType::Forward);
        self.stencils.apply(&mut elec_x, stencil::StencilType::GradXDelSquaredInv, stencil::DifferenceType::Forward);

        self.elec_y = Some(elec_y);
        self.elec_x = Some(elec_x);
        self.mag_z = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));
        self.mag_z_integral = Some(Array2::<f64>::zeros((self.geometry.nx, self.geometry.ny)));

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
        // Compute the current density from the charge density and the charge density at the last tick.
        // 
        // This is calculated by finding an approximate straight line between the locations of the charges
        
        let jx = self.current_x.as_mut().unwrap();
        let jy = self.current_y.as_mut().unwrap();
        
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
        let elec_x = self.elec_x.as_mut().unwrap();
        let elec_y = self.elec_y.as_mut().unwrap();
        let mag_z = self.mag_z.as_mut().unwrap();
        let mag_z_integral = self.mag_z_integral.as_mut().unwrap();
        let current_x = self.current_x.as_ref().unwrap();
        let current_y = self.current_y.as_ref().unwrap();


        // Evolve Bz field first, to get it from half a tick behind to half a tick ahead of the E field
        for (i, j, sigma_x, sigma_y) in pml::pml_iterator_from_geometry(&self.geometry) {
            
            let d_elec_x_dy = self.stencils.evaluate(&elec_x, i, j, &stencil::StencilType::GradY, &stencil::DifferenceType::Forward);
            let d_elec_y_dx = self.stencils.evaluate(&elec_y, i, j, &stencil::StencilType::GradX, &stencil::DifferenceType::Forward);

            // dB_z/dt = -dE_x/dy + dE_y/dx
            mag_z[[i,j]] += (d_elec_x_dy - d_elec_y_dx) * delta_t;
            
            // Now add the PML damping terms. This is a non-physical term that damps the fields near the boundary to 
            // mimic vacuum BCs. This is explained in the following sources:
            //  * https://onlinelibrary.wiley.com/doi/book/10.1002/9781118646700, chapter 3 (not very clearly written, but
            //    the equations are there)
            //  * https://arxiv.org/abs/2108.05348 (this is nicely clear but not explicit for the EM case)

            mag_z[[i,j]] += (-(sigma_x+sigma_y) * mag_z[[i, j]] - sigma_x*sigma_y * mag_z_integral[[i,j]]) * delta_t;

            // Finally update the integral of Bz over time (which is just used for the PML damping, not anything physical)
            // Note that in tests I found neglecting the integral terms in the PML didn't make a huge difference to the
            // level of reflections (just a qualitative observation, not a rigorous test). Without the integral terms,
            // it's pretty obviously just a local damping term.
            
            mag_z_integral[[i,j]] += mag_z[[i,j]] * delta_t;
            

            
        }

        // Now the B field is half a tick ahead, so update the E field a tick to get ahead again
        for (i, j, sigma_x, sigma_y) in pml::pml_iterator_from_geometry(&self.geometry) {
            let d_magz_dy = self.stencils.evaluate(&mag_z, i, j, &stencil::StencilType::GradY, &stencil::DifferenceType::Backward);
            let d_magz_dx = self.stencils.evaluate(&mag_z, i, j, &stencil::StencilType::GradX, &stencil::DifferenceType::Backward);
            elec_x[[i,j]] += (d_magz_dy - current_x[[i,j]]) * delta_t;
            elec_y[[i,j]] += (-d_magz_dx - current_y[[i,j]]) * delta_t;
            
            // PML damping, as above:
            let pml_x_term = -sigma_y * elec_x[[i, j]] + sigma_x * self.stencils.evaluate(&mag_z_integral, i, j, &stencil::StencilType::GradY, &stencil::DifferenceType::Backward);
            let pml_y_term = -sigma_x * elec_y[[i, j]] - sigma_y * self.stencils.evaluate(&mag_z_integral, i, j, &stencil::StencilType::GradX, &stencil::DifferenceType::Backward) ;
            elec_x[[i,j]] += pml_x_term*delta_t;
            elec_y[[i,j]] += pml_y_term*delta_t;
        }

    }
}

impl FieldConfiguration {
    pub fn evaluate_elec_interpolated(&mut self, x: f64, y: f64) -> (f64, f64) {
        self.ensure_initialized();
        let elec_x = evaluate_grid_interpolated_or_0(self, &self.elec_x, x, y);
        let elec_y = evaluate_grid_interpolated_or_0(self, &self.elec_y, x, y);
        (elec_x, elec_y)
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
    let mut potential: f64 = 0.0;
    for charge in &field_configuration.charges {
        let dx = x - charge.x;
        let dy = y - charge.y;
        let r = (dx * dx + dy * dy + SOFTEN).sqrt();
        potential -= FIELD_SCALING * charge.charge / r;
    }
    potential
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
    // Here we take advantage of the fact the Biot-Savart law in a plane looks like $\hat{z} \times \nabla \phi$
    // where $\phi$ is a fictional electrostatic potential generated by replacing currents with charges.
    let Pair { u, v } = compute_field_electrostatic_direct(field_configuration, x, y);
    Pair { u: v, v: -u }
}

#[wasm_bindgen]
pub fn compute_electric_field_dynamic(field_configuration: &mut FieldConfiguration, x: f64, y: f64) -> Pair {
    let (u, v) = field_configuration.evaluate_elec_interpolated(x, y);
    Pair { u, v }
}
