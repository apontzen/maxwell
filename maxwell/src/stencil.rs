use ndarray::{Array2};
use num_complex::Complex;
use crate::fourier;

pub struct FourierStencils {
    pub del_squared: Array2<Complex<f64>>,
    pub grad_x: Array2<Complex<f64>>,
    pub grad_y: Array2<Complex<f64>>,
    pub grad_x_del_squared_inv: Array2<Complex<f64>>,
    pub grad_y_del_squared_inv: Array2<Complex<f64>>,
}

pub enum StencilType {
    DelSquared,
    GradX,
    GradY,
    GradXDelSquaredInv,
    GradYDelSquaredInv,
}

impl FourierStencils {
    pub fn new(size_x: f64, size_y: f64, nx: usize, ny: usize) -> FourierStencils {
        let mut del_squared = make_del_squared_stencil(size_x, size_y, nx, ny);
        fourier::array_fft(&mut del_squared);

        let mut grad_x = make_grad_x_stencil(size_x, size_y, nx, ny);
        fourier::array_fft(&mut grad_x);
        let mut grad_y = make_grad_y_stencil(size_x, size_y, nx, ny);
        fourier::array_fft(&mut grad_y);

        fourier::array_fft_renormalise(&mut grad_x);
        fourier::array_fft_renormalise(&mut grad_y);
        
        let mut del_squared_inv = del_squared.clone();
        invert_fourier_stencil(&mut del_squared);
        
        // let mut del_squared_inv = make_1_over_r_stencil(size_x, size_y, nx, ny);
        // fourier::array_fft(&mut del_squared_inv);

        let grad_x_del_squared_inv = (&grad_x) * (&del_squared_inv);
        let grad_y_del_squared_inv = (&grad_y) * (&del_squared_inv);
        
        FourierStencils { del_squared, grad_x, grad_y, grad_x_del_squared_inv, grad_y_del_squared_inv }
    }

    pub fn apply(&self, array: &mut Array2<Complex<f64>>, stencil_type: StencilType) {
        fourier::array_fft(array);
        match stencil_type {
            StencilType::DelSquared => (*array) *= &self.del_squared,
            StencilType::GradX => (*array) *= &self.grad_x,
            StencilType::GradY => (*array) *= &self.grad_y,
            StencilType::GradXDelSquaredInv => (*array) *= &self.grad_x_del_squared_inv,
            StencilType::GradYDelSquaredInv => (*array) *= &self.grad_y_del_squared_inv,
        }
        (*array) *= Complex::new(1.0, 0.0);
        fourier::array_ifft(array);
    }

}

fn invert_fourier_stencil(array: &mut ndarray::Array2<Complex<f64>>) {
    let mut slice = array.as_slice_mut().unwrap();
    let max_abs_val = slice.iter().map(|x| x.norm()).max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)).unwrap_or(0.0);
    if max_abs_val == 0.0 {
        return;
    } else {
        slice.iter_mut().for_each(|x| *x = 
            if x.norm()>max_abs_val*1e-8 {1.0/(*x)} else {Complex::new(0.0, 0.0)});
    }
}

pub fn make_circle_stencil(radius: f64, size_x: f64, size_y: f64, nx: usize, ny: usize) -> Array2<Complex<f64>> {
    let mut E = Array2::<Complex<f64>>::zeros((nx, ny));
    let dx = size_x / nx as f64;
    let dy = size_y / ny as f64;
    let pi = std::f64::consts::PI;
    let norm = 1.0/(pi*radius*radius);

    let nx = nx as i32;
    let ny = ny as i32;
    
    for i in 0..nx {
        for j in 0..ny { 
            let i_zero_centred = (i+nx/2)%nx - nx/2;
            let j_zero_centred = (j+ny/2)%ny - ny/2;
            let x = (i_zero_centred as f64) * dx + dx/2.0;
            let y = (j_zero_centred as f64) * dy + dy/2.0;
            if x * x + y * y < radius * radius {
                E[[i as usize, j as usize]] = 50000.0* Complex::new(norm, 0.0);
            }
        }
    }
    E
}

pub fn make_1_over_r_stencil(size_x: f64, size_y: f64, nx: usize, ny: usize) -> Array2<Complex<f64>> {
    let mut E = Array2::<Complex<f64>>::zeros((nx, ny));
    let dx = size_x / nx as f64;
    let dy = size_y / ny as f64;
    let pi = std::f64::consts::PI;
    let norm = 1.0;

    let nx = nx as i32;
    let ny = ny as i32;
    
    for i in 0..nx {
        for j in 0..ny { 
            let i_zero_centred = (i+nx/2)%nx - nx/2;
            let j_zero_centred = (j+ny/2)%ny - ny/2;
            let x = (i_zero_centred as f64) * dx + dx/2.0;
            let y = (j_zero_centred as f64) * dy + dy/2.0;
            let r = (x*x + y*y).sqrt();
            if r > 0.0 {
                E[[i as usize, j as usize]] = Complex::new(norm/r, 0.0);
            }
        }
    }
    E
}

pub fn make_del_squared_stencil(size_x: f64, size_y: f64, nx: usize, ny: usize) -> Array2<Complex<f64>> {
    let mut result = Array2::<Complex<f64>>::zeros((nx, ny));
    let dx = size_x / nx as f64;
    let dy = size_y / ny as f64;

    let dx2_inv = Complex::new(1./(dx*dx), 0.0);
    let dy2_inv = Complex::new(1./(dy*dy), 0.0);

    
    result[[0,0]] = -2.0*dx2_inv + -2.0*dy2_inv;
    result[[nx-1,0]] = dx2_inv;
    result[[1,0]] = dx2_inv;
    result[[0,ny-1]] = dy2_inv;
    result[[0,1]] = dy2_inv;
    
    result
    
}

pub fn make_grad_x_stencil(size_x: f64, size_y: f64, nx: usize, ny: usize) -> Array2<Complex<f64>> {
    let mut result = Array2::<Complex<f64>>::zeros((nx, ny));
    let dx = size_x / nx as f64;
    let dy = size_y / ny as f64;

    let dx_inv_by_2 = Complex::new(0.5/(dx), 0.0);
    
    result[[1,0]] = dx_inv_by_2;
    result[[nx-1,0]] = -dx_inv_by_2;

    result
}

pub fn make_grad_y_stencil(size_x: f64, size_y: f64, nx: usize, ny: usize) -> Array2<Complex<f64>> {
    let mut result = Array2::<Complex<f64>>::zeros((nx, ny));
    let dx = size_x / nx as f64;
    let dy = size_y / ny as f64;

    let dy_inv_by_2 = Complex::new(0.5/(dy), 0.0);
    
    result[[0,1]] = dy_inv_by_2;
    result[[0,ny-1]] = -dy_inv_by_2;

    result
}


