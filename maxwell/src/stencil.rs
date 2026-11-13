use ndarray::{Array2};
use num_complex::Complex;
use crate::fourier;
use crate::geometry::Geometry;

/// The smoothing Gaussian applied to charges and currents has sigma = 2 cells, truncated at 4 sigma.
pub const SMOOTHING_SIGMA_CELLS: f64 = 2.0;
const SMOOTHING_HALF_WIDTH: isize = 8;

pub struct Stencils {
    pub del_squared_inv: Option<Array2<Complex<f64>>>,
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
        let mut s = Stencils { del_squared_inv: None, geometry };
        s.init_inv_laplacian();
        s
    }

    fn init_inv_laplacian(&mut self) {
        let mut del_squared = self.make_laplacian_stencil();
        fourier::array_fft(&mut del_squared);
        let mut del_squared_inv = del_squared.clone();
        self.invert_fourier_stencil(&mut del_squared_inv);
        fourier::array_fft_renormalise(&mut del_squared_inv);
        self.del_squared_inv = Some(del_squared_inv);
    }

    /// The gradient (forward differences) of del^-2 rho: (d/dx, d/dy) of the periodic solution of
    /// del^2 phi = rho with the k = 0 mode dropped. One FFT pair for both components.
    pub fn grad_del_squared_inv(&self, rho: &Array2<f64>) -> (Array2<f64>, Array2<f64>) {
        let mut phi = rho.clone();
        self.apply_fourier_stencil(&mut phi, self.del_squared_inv.as_ref().unwrap());
        let mut gx = Array2::<f64>::zeros(rho.dim());
        let mut gy = Array2::<f64>::zeros(rho.dim());
        self.apply_grad(&phi, &mut gx, StencilType::GradX, DifferenceType::Forward);
        self.apply_grad(&phi, &mut gy, StencilType::GradY, DifferenceType::Forward);
        (gx, gy)
    }

    fn apply_fourier_stencil(&self, array: &mut Array2<f64>, fourier_stencil: &Array2<Complex<f64>>) {
        let mut result = Array2::<Complex<f64>>::zeros(array.dim());
        result.zip_mut_with(array, |r, &a| *r = Complex::new(a, 0.0));
        fourier::array_fft(&mut result);
        (result) *= fourier_stencil;
        fourier::array_ifft(&mut result);
        array.zip_mut_with(&result, |a, &r| *a = r.re);
    }


    fn invert_fourier_stencil(&self, array: &mut ndarray::Array2<Complex<f64>>) {
        let slice = array.as_slice_mut().unwrap();
        let max_abs_val = slice.iter().map(|x| x.norm()).max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)).unwrap_or(0.0);
        if max_abs_val == 0.0 {
            return;
        } else {
            slice.iter_mut().for_each(|x| *x = 
                if x.norm()>max_abs_val*1e-8 {1.0/(*x)} else {Complex::new(0.0, 0.0)});
        }
    }

    pub fn make_laplacian_stencil(&self) -> Array2<Complex<f64>> {
        let mut result = Array2::<Complex<f64>>::zeros((self.geometry.nx, self.geometry.ny));
        let dx = self.geometry.delta_x();
        let dy = self.geometry.delta_y();

        let dx2_inv = Complex::new(1./(dx*dx), 0.0);
        let dy2_inv = Complex::new(1./(dy*dy), 0.0);

        
        result[[0,0]] = -2.0*dx2_inv + -2.0*dy2_inv;
        result[[self.geometry.nx-1,0]] = dx2_inv;
        result[[1,0]] = dx2_inv;
        result[[0,self.geometry.ny-1]] = dy2_inv;
        result[[0,1]] = dy2_inv;

        result
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
