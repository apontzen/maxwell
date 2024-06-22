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
        let i = (self.x / geometry.x_max * geometry.nx as f64) as usize;
        let j = (self.y / geometry.y_max * geometry.ny as f64) as usize;
        (i, j)
    }
}


impl Debug for Charge {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "Charge {{ x: {}, y: {}, charge: {} }}", self.x, self.y, self.charge)
    }
}


pub struct Geometry {
    x_max: f64,
    y_max: f64,
    nx: usize,
    ny: usize,
}

impl Geometry {
    fn delta_x(&self) -> f64 {
        self.x_max / self.nx as f64
    }

    fn delta_y(&self) -> f64 {
        self.y_max / self.ny as f64
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

pub fn evaluate_grid(field: &Array2<f64>, x: i32, y: i32) -> f64 {
    // Evaluate the field at the specified grid cell, or return 0 if x>=nx, y>=ny, x<0 or y<0
    let nx = field.shape()[0] as i32;
    let ny = field.shape()[1] as i32;
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

    let dx = field_config.geometry.x_max / nx as f64;
    let dy = field_config.geometry.y_max / ny as f64;

    assert!(nx == field_config.geometry.nx && ny == field_config.geometry.ny, "Grid size mismatch");

    if x < 0.0 || x > field_config.geometry.x_max || y < 0.0 || y > field_config.geometry.y_max {
        console::log_2(
            &"Point out of bounds:".into(),
            &format!("({}, {})", x, y).into()
        );
        return 0.0;
    }

    let i0: i32 = (x / field_config.geometry.x_max * nx as f64) as i32;
    let j0: i32 = (y / field_config.geometry.y_max * ny as f64) as i32;

    // work out centre of the cell we just landed in
    let x0 = (i0 as f64 + 0.5) * field_config.geometry.x_max / nx as f64;
    let y0 = (j0 as f64 + 0.5) * field_config.geometry.y_max / ny as f64;

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
        let cell_area = x_max * y_max / (nx * ny) as f64;
        let charge_normalization = 4000.0 / cell_area;
        FieldConfiguration { charges: vec![], charges_at_last_tick: vec![], 
            geometry: Geometry{x_max: x_max, y_max: y_max, nx: nx, ny: ny}, 
            cic_grid: None, Ex_grid: None, Ey_grid: None, Bz_grid: None, Bz_integral_grid: None, jx_grid: None, jy_grid: None,
            stencils: stencil::Stencils::new(x_max, y_max, nx, ny), charge_normalization }
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
                grid[[i, j]] += charge.charge * self.charge_normalization; // CIC charge 
            }

            self.stencils.apply(grid, stencil::StencilType::Soften, stencil::DifferenceType::Central);
        }
    }

    pub fn initialize_on_constraints(&mut self) {
        self.ensure_cic_grid();

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

    pub fn ensure_cic_grid(&mut self) {
        if self.cic_grid.is_none() {
            self.make_cic_grid();
        }
    }

    pub fn ensure_initialized(&mut self) {
        if self.Ey_grid.is_none() {
            self.initialize_on_constraints();
        }
    }

    pub fn evaluate_cic_grid_interpolated(&mut self, x: f64, y: f64) -> f64 {
        self.ensure_cic_grid();
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
                jx[[i, j_earlier]] = x_current_density;
            }

            for j in j_min..j_max {
                jy[[i_now, j]] = y_current_density;
            }
            
            // Now apply the stencils to smooth the current density field, in just the same way as the charge density field is smoothed
            self.stencils.apply(jx, stencil::StencilType::Soften, stencil::DifferenceType::Central);
            self.stencils.apply(jy, stencil::StencilType::Soften, stencil::DifferenceType::Central);


        }
        

    }
    
    pub fn tick(&mut self, delta_t: f64) {
        self.ensure_initialized();
        self.make_currents(delta_t);
        self.charges_at_last_tick = self.charges.clone();

        // Evolve the fields by one timestep. Note that the B and density fields are half a tick behind and half a grid cell
        // to the left of the E and j fields.

        // Evolve Bz field first.
        let Ex = self.Ex_grid.as_mut().unwrap();
        let Ey = self.Ey_grid.as_mut().unwrap();
        let Bz = self.Bz_grid.as_mut().unwrap();
        let Bz_integral = self.Bz_integral_grid.as_mut().unwrap();
        let jx = self.jx_grid.as_ref().unwrap();
        let jy = self.jy_grid.as_ref().unwrap();


        // add PML damping to the Bz field
        for (i, j, sigma_x, sigma_y) in pml::pml_iterator_from_geometry(&self.geometry) {
            Bz[[i,j]] += (-(sigma_x+sigma_y) * Bz[[i, j]])*delta_t; // + sigma_x*sigma_y * Bz_integral[[i,j]]) * delta_t;
        }


        let mut d_Ex_dy: Array2<f64> = self.stencils.apply_non_destructively(&Ex, stencil::StencilType::GradY, stencil::DifferenceType::Forward);
        let mut d_Ey_dx: Array2<f64> = self.stencils.apply_non_destructively(&Ey, stencil::StencilType::GradX, stencil::DifferenceType::Forward);

        d_Ex_dy *= delta_t;
        d_Ey_dx *= delta_t;

        (*Bz) += &d_Ex_dy;
        (*Bz) -= &d_Ey_dx;

        
        *Bz_integral += &(Bz as &_ * delta_t);

        


        // evolve Ex and Ey fields
        let mut d_Bz_dy = self.stencils.apply_non_destructively(&Bz, stencil::StencilType::GradY, stencil::DifferenceType::Backward);
        let mut d_Bz_dx = self.stencils.apply_non_destructively(&Bz, stencil::StencilType::GradX, stencil::DifferenceType::Backward);


        d_Bz_dy *= delta_t;
        d_Bz_dx *= delta_t;

        (*Ex) += &d_Bz_dy;
        
        (*Ey) -= &d_Bz_dx;

        (*Ex) -= &(jx*delta_t);
        (*Ey) -= &(jy*delta_t);

        // add PML damping to the Ex and Ey fields
        for (i, j, sigma_x, sigma_y) in pml::pml_iterator_from_geometry(&self.geometry) {
            let pml_x_term = -sigma_y * Ex[[i, j]]; // + sigma_x * self.stencils.evaluate(&Bz_integral, i, j, &stencil::StencilType::GradY, &stencil::DifferenceType::Backward);
            let pml_y_term = -sigma_x * Ey[[i, j]]; //- sigma_y * self.stencils.evaluate(&Bz_integral, i, j, &stencil::StencilType::GradX, &stencil::DifferenceType::Backward) ;
            Ex[[i,j]] += pml_x_term*delta_t;
            Ey[[i,j]] += pml_y_term*delta_t;
        }

        
    
        

    }
}

impl FieldConfiguration {
    // For internal functions only

    pub fn evaluate_E_grid_interpolated(&mut self, x: f64, y: f64) -> (f64, f64) {
        self.ensure_initialized();
        let Ex = evaluate_grid_interpolated_or_0(self, &self.Ex_grid, x, y);
        let Ey = evaluate_grid_interpolated_or_0(self, &self.Ey_grid, x, y);
        (Ex, Ey)
    }

    pub fn count_nans(name: &str, grid: &Array2<f64>) {
        let mut nan_count = 0;
        for element in grid.iter() {
            if element.is_nan() || element.is_infinite() {
                nan_count += 1;
            }
        }
        web_sys::console::log_1(&format!("Number of NaN values in {}: {}", name, nan_count).into());
    }

}


#[wasm_bindgen]
pub fn compute_field_electrostatic_direct(field_configuration: &FieldConfiguration, x: f64, y: f64) -> JsValue {
    let mut u: f64 = 0.0;
    let mut v: f64 = 0.0;
    for charge in &field_configuration.charges {
        let dx = x - charge.x;
        let dy = y - charge.y;
        let r = (dx * dx + dy * dy).sqrt();
        let k = 20000.0 * charge.charge;
        u += k * dx / (r * r * r);
        v += k * dy / (r * r * r);
    }
    to_value(&(u, v)).unwrap()
}

#[wasm_bindgen]
pub fn compute_field_magnetostatic_direct(field_configuration: &FieldConfiguration, x: f64, y:f64) -> JsValue {
    let mut u: f64 = 0.0;
    let mut v: f64 = 0.0;
    for current in &field_configuration.charges { // charges interpreted as current
        let dx = x - current.x;
        let dy = y - current.y;
        let r = (dx * dx + dy * dy).sqrt();
        let k = 20000.0 * current.charge;
        u += k * dy / (r * r * r);
        v -= k * dx / (r * r * r);
    }
    to_value(&(u, v)).unwrap()
}

#[wasm_bindgen]
pub fn compute_field_electrostatic_fourier(field_configuration: &mut FieldConfiguration, x: f64, y: f64) -> JsValue {
    let (u, v) = field_configuration.evaluate_E_grid_interpolated(x, y);
    to_value(&(u, v)).unwrap()
}