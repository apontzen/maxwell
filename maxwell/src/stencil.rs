use ndarray::{Array2};
use num_complex::Complex;
use crate::fourier;

pub struct FourierStencils {
    pub del_squared_inv: Array2<Complex<f64>>,
   
    size_x: f64,
    size_y: f64,
    nx: usize,
    ny: usize,

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
        let mut del_squared_inv = del_squared.clone();
        invert_fourier_stencil(&mut del_squared_inv);
        fourier::array_fft_renormalise(&mut del_squared_inv);

        FourierStencils { del_squared_inv, size_x, size_y, nx, ny}
    }

    pub fn apply_inv_laplacian(&self, array: &mut Array2<Complex<f64>>) -> Array2<Complex<f64>> {
        let mut result = array.clone();
        fourier::array_fft(&mut result);
        (result) *= &self.del_squared_inv;
        fourier::array_ifft(&mut result);
        result
    }

    pub fn apply(&self, array: &mut Array2<Complex<f64>>, stencil_type: StencilType) {
        
        let mut array_temp = match stencil_type {
            StencilType::GradXDelSquaredInv | StencilType::GradYDelSquaredInv => self.apply_inv_laplacian(array),
            _ => array.clone()
        };
        
        
        match stencil_type {
            StencilType::GradXDelSquaredInv | StencilType::GradX => array.clone_from(&apply_grad_x_stencil(&array_temp, self.size_x/self.nx as f64)),
            StencilType::GradYDelSquaredInv | StencilType::GradY => array.clone_from(&apply_grad_y_stencil(&array_temp, self.size_y/self.ny as f64)),
            _ => array.clone_from(&array_temp),
        };
        
        
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

pub fn apply_grad_x_stencil(array: & Array2<Complex<f64>>, dx: f64) -> Array2<Complex<f64>> {
    let mut result = Array2::<Complex<f64>>::zeros(array.dim());
    let dx_inv_by_2 = Complex::new(0.5/(dx), 0.0);
    let nx = array.dim().0;
    let ny = array.dim().1;
    for i in 0..nx {
        for j in 0..ny {
            let i_plus_one = if(i<nx-1) {i+1} else {0};
            let i_minus_one = if(i>0) {i-1} else {nx-1};
            result[[i,j]] = dx_inv_by_2 * (-array[[i_plus_one, j]] + array[[i_minus_one, j]]);
        }
    }
    result
}

pub fn apply_grad_y_stencil(array: & Array2<Complex<f64>>, dy: f64) -> Array2<Complex<f64>> {
    let mut result = Array2::<Complex<f64>>::zeros(array.dim());
    let dy_inv_by_2 = Complex::new(0.5/(dy), 0.0);
    let nx = array.dim().0;
    let ny = array.dim().1;
    for i in 0..nx {
        for j in 0..ny {
            let j_plus_one = if j < ny - 1 { j + 1 } else { 0 };
            let j_minus_one = if j > 0 { j - 1 } else { ny - 1 };
            result[[i, j]] = dy_inv_by_2 * (-array[[i, j_plus_one]] + array[[i, j_minus_one]]);
        }
    }
    result
}