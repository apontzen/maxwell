// 2D FFTs of ndarray arrays (unnormalised, as rustfft), with the FFT plans cached per size: planning
// is expensive for sizes with large prime factors (the default grid is 111 = 3 x 37) and these
// transforms run every frame.

use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::Arc;
use ndarray::Array2;
use num_complex::Complex;
use rustfft::{Fft, FftPlanner};

struct Plans {
    planner: FftPlanner<f64>,
    cache: HashMap<(usize, bool), Arc<dyn Fft<f64>>>,
}

thread_local! {
    static PLANS: RefCell<Plans> = RefCell::new(Plans { planner: FftPlanner::new(), cache: HashMap::new() });
}

fn plan(n: usize, inverse: bool) -> Arc<dyn Fft<f64>> {
    PLANS.with(|p| {
        let mut p = p.borrow_mut();
        if let Some(f) = p.cache.get(&(n, inverse)) {
            return f.clone();
        }
        let f = if inverse { p.planner.plan_fft_inverse(n) } else { p.planner.plan_fft_forward(n) };
        p.cache.insert((n, inverse), f.clone());
        f
    })
}

fn fft_2d(array: &mut Array2<Complex<f64>>, inverse: bool) {
    let (nx, ny) = array.dim();
    let rows = plan(ny, inverse);
    let cols = plan(nx, inverse);
    let mut scratch = vec![Complex::new(0.0, 0.0); rows.get_inplace_scratch_len().max(cols.get_inplace_scratch_len())];
    {
        // the rows of a standard-layout array are contiguous
        let data = array.as_slice_mut().expect("standard layout");
        for row in data.chunks_exact_mut(ny) {
            rows.process_with_scratch(row, &mut scratch);
        }
    }
    let mut column = vec![Complex::new(0.0, 0.0); nx];
    for j in 0..ny {
        for i in 0..nx { column[i] = array[[i, j]]; }
        cols.process_with_scratch(&mut column, &mut scratch);
        for i in 0..nx { array[[i, j]] = column[i]; }
    }
}

pub fn array_fft(array: &mut Array2<Complex<f64>>) {
    fft_2d(array, false);
}

pub fn array_ifft(array: &mut Array2<Complex<f64>>) {
    fft_2d(array, true);
}

pub fn array_fft_renormalise(array: &mut Array2<Complex<f64>>) {
    let norm_factor = (array.len()) as f64;
    array.mapv_inplace(|x| x / norm_factor);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_the_definition_and_inverts() {
        let (nx, ny) = (6, 5);
        let a = Array2::from_shape_fn((nx, ny), |(i, j)| Complex::new((i * 7 + j * 3) as f64 % 5.0 - 2.0, (i + 2 * j) as f64 % 3.0));
        let mut b = a.clone();
        array_fft(&mut b);
        for k in 0..nx {
            for l in 0..ny {
                let mut s = Complex::new(0.0, 0.0);
                for i in 0..nx {
                    for j in 0..ny {
                        let ph = -2.0 * std::f64::consts::PI * ((k * i) as f64 / nx as f64 + (l * j) as f64 / ny as f64);
                        s += a[[i, j]] * Complex::from_polar(1.0, ph);
                    }
                }
                assert!((s - b[[k, l]]).norm() < 1e-9);
            }
        }
        array_ifft(&mut b);
        array_fft_renormalise(&mut b);
        for (x, y) in a.iter().zip(b.iter()) { assert!((x - y).norm() < 1e-12); }
    }
}
