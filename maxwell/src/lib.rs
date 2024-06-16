use wasm_bindgen::prelude::*;
use serde::{Serialize, Deserialize};
use serde_wasm_bindgen::{from_value, to_value};
use console_error_panic_hook;
use web_sys::console as console;
use ndarray::Array2;
use std::fmt::{self, Debug, Formatter};

use fft2d::{fft_2d, fftshift};

#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
    console::log_1(&"Hello from Rust, panic hook has been initialised!".into());
}

#[wasm_bindgen]
#[derive(Serialize, Deserialize)]
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

impl Debug for Charge {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "Charge {{ x: {}, y: {}, charge: {} }}", self.x, self.y, self.charge)
    }
}

#[wasm_bindgen]
pub struct FieldConfiguration {
    charges: Vec<Charge>,
    x_max: f64,
    y_max: f64,
    nx: usize,
    ny: usize,
    cic_grid: Option<Array2<Complex<f64>>>,
    E_grid: Option<Array2<Complex<f64>>>,
}

pub fn evaluate_grid(field: &Array2<Complex<f64>>, x: i32, y: i32) -> f64 {
    // Evaluate the field at the specified grid cell, or return 0 if x>=nx, y>=ny, x<0 or y<0
    let nx = field.shape()[0] as i32;
    let ny = field.shape()[1] as i32;
    if x < 0 || x >= nx || y < 0 || y >= ny {
        return 0.0;
    }
    field[[x as usize, y as usize]].re
}

pub fn evaluate_grid_interpolated_or_0(field_config: &FieldConfiguration, grid: &Option<Array2<Complex<f64>>>, x: f64, y: f64) -> f64 {
    if let Some(grid) = grid {
        evaluate_grid_interpolated(field_config, grid, x, y)
    } else {
        0.0
    }
}

pub fn evaluate_grid_interpolated(field_config: &FieldConfiguration, grid: &Array2<Complex<f64>>, x: f64, y: f64) -> f64 {
    let nx = grid.shape()[0];
    let ny = grid.shape()[1];

    let dx = field_config.x_max / nx as f64;
    let dy = field_config.y_max / ny as f64;

    assert!(nx == field_config.nx && ny == field_config.ny, "Grid size mismatch");

    if x < 0.0 || x > field_config.x_max || y < 0.0 || y > field_config.y_max {
        console::log_2(
            &"Point out of bounds:".into(),
            &format!("({}, {})", x, y).into()
        );
        return 0.0;
    }

    let i0: i32 = (x / field_config.x_max * nx as f64) as i32;
    let j0: i32 = (y / field_config.y_max * ny as f64) as i32;

    // work out centre of the cell we just landed in
    let x0 = (i0 as f64 + 0.5) * field_config.x_max / nx as f64;
    let y0 = (j0 as f64 + 0.5) * field_config.y_max / ny as f64;

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
        FieldConfiguration { charges: vec![], x_max: x_max, y_max: y_max, nx: nx, ny: ny, cic_grid: None, E_grid: None,}
    }

    pub fn set_charges(&mut self, charges: JsValue) {
        self.charges = match from_value(charges) {
            Ok(charges) => charges,
            Err(err) => {
                console::log_1(&format!("Error deserializing charges: {:?}", err).into());
                vec![]
            }
        };
        self.cic_grid = None;
    }

    pub fn make_cic_grid(&mut self) {
        self.cic_grid = Some(Array2::<Complex<f64>>::zeros((self.nx, self.ny)));
        if let Some(ref mut grid) = self.cic_grid {
            for charge in &self.charges {
                if(charge.x<0.0 || charge.x>self.x_max || charge.y<0.0 || charge.y>self.y_max) {
                    console::log_1(&format!("Charge out of bounds: {:?}", charge).into());
                    continue;
                }
                let i = (charge.x / self.x_max * self.nx as f64) as usize;
                let j = (charge.y / self.y_max * self.ny as f64) as usize;
                grid[[i, j]] += charge.charge * 20.0;

            }
        }
    }

    pub fn make_E_grid(&mut self) {
        self.ensure_cic_grid();

        let cic_grid: &mut Array2<Complex<f64>> = self.cic_grid.as_mut().unwrap();

        cic_grid.as_slice_mut()

        let E = cic_grid.as_slice_mut().unwrap().fft();

        for i in 0..self.nx {
            for j in 0..self.ny {
                let kx = if i < self.nx / 2 { i as f64 } else { (i as f64) - (self.nx as f64) };
                let ky = if j < self.ny / 2 { j as f64 } else { (j as f64) - (self.ny as f64) };
                let k = (kx * kx + ky * ky).sqrt();
                if k == 0.0 {
                    E[[i, j]] = Complex::zero();
                } else {
                    E[[i, j]] = E[[i, j]] * Complex::new(0.0, -k);
                }
            }
        }

        self.E_grid = Some(E);




    }

    pub fn ensure_cic_grid(&mut self) {
        if self.cic_grid.is_none() {
            self.make_cic_grid();
        }
    }

    pub fn ensure_E_grid(&mut self) {
        if self.E_grid.is_none() {
            self.make_E_grid();
        }
    }

    pub fn evaluate_cic_grid_interpolated(&self, x: f64, y: f64) -> f64 {
        evaluate_grid_interpolated_or_0(self, &self.cic_grid, x, y)
    }

    pub fn evaluate_E_grid_interpolated(&self, x: f64, y: f64) -> f64 {
        evaluate_grid_interpolated_or_0(self, &self.E_grid, x, y)
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
pub fn compute_field_electrostatic_fourier(field_configuration: &mut FieldConfiguration, x: f64, y: f64) -> JsValue {
    field_configuration.ensure_cic_grid();
    let u: f64 = 0.0;
    let v: f64 = field_configuration.evaluate_E_grid_interpolated(x, y);
    to_value(&(u, v)).unwrap()
}