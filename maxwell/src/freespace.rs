// The static potential of charges in free space (vacuum boundary conditions) on the padded grid.
//
// The grid solvers used to solve Poisson's equation periodically (FFT), so that the static field
// in the canvas included that of the periodic images of the charges. That distorted the field
// noticeably (the canvas is only twice as high as the separation of a typical dipole), and the
// dynamic solver relaxed to that distorted field too, since its absorbing layer pulls the field
// towards the static one.
//
// Here instead: the outermost ring of cells of the padded grid is at least the absorbing layer's
// thickness (20 cells) away from any charge, where the potential of the (Gaussian-smeared)
// charges equals that of point charges, sum Q ln(r) / 2 pi, to within ~exp(-r^2 / 2 sigma^2). So
// the potential is set to that on the ring, and the interior is solved exactly with the grid's
// 5-point Laplacian (the one whose forward-difference gradient gives E, as for the periodic solve),
// by sine transforms (DST-I). The grid's Gauss's law then holds exactly everywhere inside the
// ring, as the dynamic solver requires.

use ndarray::Array2;
use num_complex::Complex;
use crate::fourier;
use crate::geometry::Geometry;

/// DST-I of the n values x[offset + stride * k], in place: X_k = sum_j x_j sin(pi (j+1)(k+1) / (n+1)),
/// via a complex FFT of length 2(n+1) of the odd extension.
fn dst1(x: &mut [f64], offset: usize, stride: usize, n: usize, buf: &mut Vec<Complex<f64>>, scratch: &mut Vec<Complex<f64>>) {
    let len = 2 * (n + 1);
    let fft = fourier::plan(len, false);
    buf.clear();
    buf.resize(len, Complex::new(0.0, 0.0));
    for j in 0..n {
        let v = x[offset + stride * j];
        buf[j + 1] = Complex::new(v, 0.0);
        buf[len - 1 - j] = Complex::new(-v, 0.0);
    }
    let need = fft.get_inplace_scratch_len();
    if scratch.len() < need { scratch.resize(need, Complex::new(0.0, 0.0)); }
    fft.process_with_scratch(buf, scratch);
    for k in 0..n { x[offset + stride * k] = -0.5 * buf[k + 1].im; }
}

/// The potential phi of the charge density `rho` (on the padded grid), with del^2 phi = rho
/// (5-point Laplacian) inside the outermost ring of cells and the free-space potential of the
/// point charges `charges` = (x, y, Q), Q the total of each charge's deposit (rho times cell
/// area), on the ring. `rho` must vanish on and next to the ring.
pub fn potential(geometry: &Geometry, rho: &Array2<f64>, charges: &[(f64, f64, f64)]) -> Array2<f64> {
    let (nx, ny, nb) = (geometry.nx, geometry.ny, geometry.nboundary as f64);
    let (dx, dy) = (geometry.delta_x(), geometry.delta_y());
    let tau = 2.0 * std::f64::consts::PI;
    let mut phi = Array2::<f64>::zeros((nx, ny));
    // the ring
    let free = |i: usize, j: usize| -> f64 {
        let (x, y) = ((i as f64 - nb + 0.5) * dx, (j as f64 - nb + 0.5) * dy);
        charges.iter().map(|&(cx, cy, q)| q * ((x - cx).powi(2) + (y - cy).powi(2)).ln() / (2.0 * tau)).sum()
    };
    for i in 0..nx { phi[[i, 0]] = free(i, 0); phi[[i, ny - 1]] = free(i, ny - 1); }
    for j in 1..ny - 1 { phi[[0, j]] = free(0, j); phi[[nx - 1, j]] = free(nx - 1, j); }
    // the interior (nx - 2) x (ny - 2): rho minus the known ring values' part of the Laplacian
    let (mx, my) = (nx - 2, ny - 2);
    let (cx, cy) = (1.0 / (dx * dx), 1.0 / (dy * dy));
    let mut b = vec![0.0; mx * my];
    for a in 0..mx {
        for c in 0..my {
            let (i, j) = (a + 1, c + 1);
            let mut v = rho[[i, j]];
            if a == 0 { v -= cx * phi[[0, j]]; }
            if a == mx - 1 { v -= cx * phi[[nx - 1, j]]; }
            if c == 0 { v -= cy * phi[[i, 0]]; }
            if c == my - 1 { v -= cy * phi[[i, ny - 1]]; }
            b[a * my + c] = v;
        }
    }
    let (mut buf, mut scratch) = (Vec::new(), Vec::new());
    for a in 0..mx { dst1(&mut b, a * my, 1, my, &mut buf, &mut scratch); }
    for c in 0..my { dst1(&mut b, c, my, mx, &mut buf, &mut scratch); }
    // eigenvalues of the 1D second differences with Dirichlet ends
    let lam = |n: usize, c: f64| -> Vec<f64> {
        (1..=n).map(|k| -4.0 * c * (std::f64::consts::PI * k as f64 / (2.0 * (n + 1) as f64)).sin().powi(2)).collect()
    };
    let (lx, ly) = (lam(mx, cx), lam(my, cy));
    // the inverse DST-I is the DST-I times 2 / (n + 1), in each direction
    let norm = 4.0 / ((mx + 1) * (my + 1)) as f64;
    for a in 0..mx { for c in 0..my { b[a * my + c] *= norm / (lx[a] + ly[c]); } }
    for c in 0..my { dst1(&mut b, c, my, mx, &mut buf, &mut scratch); }
    for a in 0..mx { dst1(&mut b, a * my, 1, my, &mut buf, &mut scratch); }
    for a in 0..mx { for c in 0..my { phi[[a + 1, c + 1]] = b[a * my + c]; } }
    phi
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dst1_matches_its_definition_and_inverts() {
        let (mut buf, mut scratch) = (Vec::new(), Vec::new());
        for n in [1usize, 2, 7, 54] {
            let x: Vec<f64> = (0..n).map(|j| ((j * 37 + 11) % 17) as f64 - 8.0).collect();
            let mut y = x.clone();
            dst1(&mut y, 0, 1, n, &mut buf, &mut scratch);
            for k in 0..n {
                let want: f64 = (0..n).map(|j| x[j] * (std::f64::consts::PI * ((j + 1) * (k + 1)) as f64 / (n + 1) as f64).sin()).sum();
                assert!((y[k] - want).abs() < 1e-9 * (1.0 + want.abs()), "n={} k={}", n, k);
            }
            dst1(&mut y, 0, 1, n, &mut buf, &mut scratch);
            for j in 0..n { assert!((y[j] * 2.0 / (n + 1) as f64 - x[j]).abs() < 1e-9, "n={} j={}", n, j); }
        }
    }

    #[test]
    fn solves_poisson_exactly_with_the_free_space_ring() {
        let g = Geometry { x_max: 668.0, y_max: 335.0, nx: 151, ny: 96, nboundary: 20 };
        let (dx, dy) = (g.delta_x(), g.delta_y());
        // two blobs of opposite charge, and net charge
        let mut rho = Array2::<f64>::zeros((g.nx, g.ny));
        let charges = [(200.0, 150.0, 3.0), (450.0, 200.0, -1.0)];
        for &(x, y, q) in &charges {
            let (u, v) = ((x / dx + 20.0) as usize, (y / dy + 20.0) as usize);
            rho[[u, v]] += q / (dx * dy);
        }
        let phi = potential(&g, &rho, &charges);
        let mut worst: f64 = 0.0;
        for i in 1..g.nx - 1 { for j in 1..g.ny - 1 {
            let lap = (phi[[i + 1, j]] - 2.0 * phi[[i, j]] + phi[[i - 1, j]]) / (dx * dx)
                + (phi[[i, j + 1]] - 2.0 * phi[[i, j]] + phi[[i, j - 1]]) / (dy * dy);
            worst = worst.max((lap - rho[[i, j]]).abs());
        } }
        let scale = 3.0 / (dx * dy);
        assert!(worst < 1e-9 * scale, "{}", worst / scale);
    }
}
