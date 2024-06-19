use ndarray::{Array2};
use fft2d::slice::{fft_2d, ifft_2d};
use num_complex::Complex;


pub fn array_fft(array: &mut ndarray::Array2<Complex<f64>> ) {
    let width = array.shape()[0];
    let height = array.shape()[1];
    let mut slice = array.as_slice_mut().unwrap();
    fft_2d( width, height, &mut slice);
}

pub fn array_ifft(array: &mut ndarray::Array2<Complex<f64>> ) {
    let width = array.shape()[0];
    let height = array.shape()[1];
    let mut slice = array.as_slice_mut().unwrap();
    ifft_2d( width, height, &mut slice);
}

pub fn array_fft_renormalise(array: &mut ndarray::Array2<Complex<f64>> ) {
    let width = array.shape()[0];
    let height = array.shape()[1];
    let norm_factor = width*height;
    let mut slice = array.as_slice_mut().unwrap();
    slice.iter_mut().for_each(|x| *x = *x / (norm_factor as f64) );
}
