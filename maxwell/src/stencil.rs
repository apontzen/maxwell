use ndarray::{Array2};
use num_complex::Complex;
use crate::fourier;

pub struct Stencils {
    pub del_squared_inv: Option<Array2<Complex<f64>>>,
   
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

pub enum DifferenceType {
    Forward,
    Backward,
    Central,
}

impl Stencils {
    pub fn new(size_x: f64, size_y: f64, nx: usize, ny: usize) -> Stencils {
        let mut s = Stencils { del_squared_inv: None, size_x, size_y, nx, ny};
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

    fn apply_inv_laplacian(&self, array: &mut Array2<f64>) {
        let mut result = Array2::<Complex<f64>>::zeros(array.dim());
        result.zip_mut_with(array, |r, &a| *r = Complex::new(a, 0.0));
        fourier::array_fft(&mut result);
        (result) *= self.del_squared_inv.as_ref().unwrap();
        fourier::array_ifft(&mut result);
        array.zip_mut_with(&result, |a, &r| *a = r.re);
    }


    fn invert_fourier_stencil(&self, array: &mut ndarray::Array2<Complex<f64>>) {
        let mut slice = array.as_slice_mut().unwrap();
        let max_abs_val = slice.iter().map(|x| x.norm()).max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)).unwrap_or(0.0);
        if max_abs_val == 0.0 {
            return;
        } else {
            slice.iter_mut().for_each(|x| *x = 
                if x.norm()>max_abs_val*1e-8 {1.0/(*x)} else {Complex::new(0.0, 0.0)});
        }
    }

    pub fn make_laplacian_stencil(&self) -> Array2<Complex<f64>> {
        let mut result = Array2::<Complex<f64>>::zeros((self.nx, self.ny));
        let dx = self.size_x / self.nx as f64;
        let dy = self.size_y / self.ny as f64;

        let dx2_inv = Complex::new(1./(dx*dx), 0.0);
        let dy2_inv = Complex::new(1./(dy*dy), 0.0);

        
        result[[0,0]] = -2.0*dx2_inv + -2.0*dy2_inv;
        result[[self.nx-1,0]] = dx2_inv;
        result[[1,0]] = dx2_inv;
        result[[0,self.ny-1]] = dy2_inv;
        result[[0,1]] = dy2_inv;
        
        result
        
    }

    pub fn apply_grad_x_stencil(&self, array: & Array2<f64>, result: &mut Array2<f64>, difference_type: DifferenceType) {
        let dx = self.size_x / self.nx as f64;
        let dx_inv = 1.0/dx;
        let dx_inv_by_2 = 0.5/dx;
        let nx = array.dim().0;
        let ny = array.dim().1;
        for i in 0..nx {
            for j in 0..ny {
                let i_plus_one = if(i<nx-1) {i+1} else {0};
                let i_minus_one = if(i>0) {i-1} else {nx-1};
                result[[i,j]] = match(difference_type) {
                    DifferenceType::Forward => dx_inv * (array[[i_plus_one, j]] - array[[i, j]]),
                    DifferenceType::Backward => dx_inv * (array[[i, j]] - array[[i_minus_one, j]]),
                    DifferenceType::Central => dx_inv_by_2 * (array[[i_plus_one, j]] - array[[i_minus_one, j]]),
                };
            }
        }
    }

    pub fn apply_grad_y_stencil(&self, array: & Array2<f64>, result: &mut Array2<f64>, difference_type: DifferenceType)  {
        let dy = self.size_y / self.ny as f64;
        let dy_inv = 1.0/dy;
        let dy_inv_by_2 = 0.5/(dy);
        let nx = array.dim().0;
        let ny = array.dim().1;
        for i in 0..nx {
            for j in 0..ny {
                let j_plus_one = if j < ny - 1 { j + 1 } else { 0 };
                let j_minus_one = if j > 0 { j - 1 } else { ny - 1 };
                result[[i,j]] = match(difference_type) {
                    DifferenceType::Forward => dy_inv * (array[[i, j_plus_one]] - array[[i, j]]),
                    DifferenceType::Backward => dy_inv * (array[[i, j]] - array[[i, j_minus_one]]),
                    DifferenceType::Central => dy_inv_by_2 * (array[[i, j_plus_one]] - array[[i, j_minus_one]]),
                };
            }
        }
    }

    pub fn apply(&self, array: &mut Array2<f64>, stencil_type: StencilType, difference_type: DifferenceType) {
        let mut scratch = array.clone();
        
        match stencil_type {
            StencilType::GradXDelSquaredInv | StencilType::GradYDelSquaredInv => self.apply_inv_laplacian(&mut scratch),
            _ => ()
        };
        
        match stencil_type {
            StencilType::GradXDelSquaredInv | StencilType::GradX => self.apply_grad_x_stencil(&scratch, array, difference_type),
            StencilType::GradYDelSquaredInv | StencilType::GradY => self.apply_grad_y_stencil(&scratch, array, difference_type),
            _ => array.clone_from(&scratch),
        };
    }

    pub fn apply_non_destructively(&self, array: &Array2<f64>, stencil_type: StencilType, difference_type: DifferenceType) -> Array2<f64> {
        let mut result = array.clone();
        self.apply(&mut result, stencil_type, difference_type);
        result
    }

}
