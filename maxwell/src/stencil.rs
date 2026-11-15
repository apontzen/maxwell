use ndarray::{Array2};
use crate::geometry::Geometry;

/// The smoothing Gaussian applied to charges and currents has sigma = 2 cells, truncated at 4 sigma.
pub const SMOOTHING_SIGMA_CELLS: f64 = 2.0;
const SMOOTHING_HALF_WIDTH: isize = 8;

pub struct Stencils {
    geometry: Geometry,
}

#[allow(dead_code)]
pub enum StencilType {
    DelSquared,
    GradX,
    GradY,
}

#[allow(dead_code)]
pub enum DifferenceType {
    Forward,
    Backward,
    Central,
}


impl Stencils {
    pub fn new(geometry: Geometry) -> Stencils {
        Stencils { geometry }
    }

    /// The forward-difference gradient of phi (as used for the static field)
    pub fn grad_forward(&self, phi: &Array2<f64>) -> (Array2<f64>, Array2<f64>) {
        let mut gx = Array2::<f64>::zeros(phi.dim());
        let mut gy = Array2::<f64>::zeros(phi.dim());
        self.apply_grad(phi, &mut gx, StencilType::GradX, DifferenceType::Forward);
        self.apply_grad(phi, &mut gy, StencilType::GradY, DifferenceType::Forward);
        (gx, gy)
    }

    /// 1D weights of the smoothing Gaussian centred on the (sub-cell) position `pos` along an axis
    /// with cells of size `delta` (sigma = 2 cells), normalised to sum to exactly 1. Returns the
    /// padded-grid index of the first weight (may be negative or beyond the grid: wrap it).
    pub fn gaussian_weights(&self, pos: f64, delta: f64) -> (isize, Vec<f64>) {
        let u = pos / delta + self.geometry.nboundary as f64; // in cells; cell i spans [i, i + 1]
        let centre = (u - 0.5).floor() as isize;
        let start = centre - SMOOTHING_HALF_WIDTH;
        let mut w: Vec<f64> = (0..(2 * SMOOTHING_HALF_WIDTH + 2) as usize).map(|k| {
            let d = (start + k as isize) as f64 + 0.5 - u;
            (-d * d / (2.0 * SMOOTHING_SIGMA_CELLS * SMOOTHING_SIGMA_CELLS)).exp()
        }).collect();
        let total: f64 = w.iter().sum();
        w.iter_mut().for_each(|v| *v /= total);
        (start, w)
    }

    /// Add value * wx[i - sx] * wy[j - sy] (periodic wrap).
    pub fn add_separable(&self, array: &mut Array2<f64>, sx: isize, wx: &[f64], sy: isize, wy: &[f64], value: f64) {
        let (nx, ny) = (self.geometry.nx as isize, self.geometry.ny as isize);
        for (a, vx) in wx.iter().enumerate() {
            if *vx == 0.0 { continue; }
            let i = (sx + a as isize).rem_euclid(nx) as usize;
            for (b, vy) in wy.iter().enumerate() {
                let j = (sy + b as isize).rem_euclid(ny) as usize;
                array[[i, j]] += value * vx * vy;
            }
        }
    }

    pub fn apply_grad(&self, array: & Array2<f64>, result: &mut Array2<f64>, stencil_type: StencilType, difference_type: DifferenceType) {
        for i in 0..self.geometry.nx {
            for j in 0..self.geometry.ny {
                result[[i,j]] = self.evaluate(array, i, j, &stencil_type, &difference_type);
            }
        }
    }

    pub fn evaluate(&self, array: &Array2<f64>, i: usize, j: usize, stencil_type: &StencilType, difference_type: &DifferenceType) -> f64 {
        let dx = self.geometry.delta_x();
        let dx_inv = 1.0/dx;
        let dx_inv_by_2 = 0.5/dx;
        let dy = self.geometry.delta_y();
        let dy_inv = 1.0/dy;
        let dy_inv_by_2 = 0.5/(dy);
        let nx = array.dim().0;
        let ny = array.dim().1;
        match stencil_type {
            StencilType::GradX => {
                let i_plus_one = if i<nx-1 {i+1} else {0};
                let i_minus_one = if i>0 {i-1} else {nx-1};
                match difference_type {
                    DifferenceType::Forward => dx_inv * (array[[i_plus_one, j]] - array[[i, j]]),
                    DifferenceType::Backward => dx_inv * (array[[i, j]] - array[[i_minus_one, j]]),
                    DifferenceType::Central => dx_inv_by_2 * (array[[i_plus_one, j]] - array[[i_minus_one, j]]),
                }
            }, 
            StencilType::GradY => {
                let j_plus_one = if j < ny - 1 { j + 1 } else { 0 };
                let j_minus_one = if j > 0 { j - 1 } else { ny - 1 };
                match difference_type {
                    DifferenceType::Forward => dy_inv * (array[[i, j_plus_one]] - array[[i, j]]),
                    DifferenceType::Backward => dy_inv * (array[[i, j]] - array[[i, j_minus_one]]),
                    DifferenceType::Central => dy_inv_by_2 * (array[[i, j_plus_one]] - array[[i, j_minus_one]]),
                }
            },
            _ => { panic!("Invalid stencil type for real-space evaluation") }
        }
    }

}
