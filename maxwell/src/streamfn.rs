// Stream function of the grid electric field, for drawing field lines of the grid solvers.
// See https://github.com/apontzen/maxwell/issues/3 (stage 2b).
//
// In 2D, field lines are contours of a stream function psi with grad psi = rot90(E) / (flux per
// unit charge). Around a charge q psi changes by q, so psi is multivalued; we split
//     psi = psi_sing + psi_rest,
// where psi_sing = sum_i q_i theta_i / 2 pi is the stream function of point charges at the
// charges' true positions (handled analytically, multivalued, in JS) and psi_rest is single valued.
// psi_rest is the least-squares fit to the remaining flux through every grid face of the physical
// region. It represents the images of the periodic solver, the curl of a time-dependent field and
// the difference between the smeared / cell-snapped charges and point charges at their true
// positions (exactly, outside the smearing kernel).
//
// The grid is staggered so that the discrete Gauss law is div_backward E = rho:
//   elec_x[i, j] is the flux density through the face between cells (i, j) and (i+1, j);
//   elec_y[i, j] through the face between cells (i, j) and (i, j+1).
// psi lives on cell corners. Corner (a, b) is the lower-left corner of cell (a, b); the corners of
// the physical region are a, b = nboundary ..= nboundary + nphys, at canvas positions
// ((a - nboundary) dx, (b - nboundary) dy).

use std::sync::Arc;
use ndarray::Array2;
use num_complex::Complex;
use rustfft::{Fft, FftPlanner};
use crate::geometry::Geometry;

/// Regularisation length of the subtracted point charges, in cells, and the distance (in units of
/// it) beyond which the regularisation is negligible
const REG_CELLS: f64 = 1.0;
const REG_FAR: f64 = 4.5;   // exp(-REG_FAR^2 / 2) ~ 4e-5
/// Within this distance (in units of the regularisation length) of a charge a segment's d theta
/// is integrated directly; beyond, the exact d theta is corrected (a segment is at most one length long)
const REG_NEAR: f64 = 2.0;

/// Gauss-Legendre nodes and weights on [-1, 1]
const GAUSS_LEGENDRE_4: [(f64, f64); 4] = [
    (-0.8611363115940526, 0.3478548451374538), (-0.3399810435848563, 0.6521451548625461),
    (0.3399810435848563, 0.6521451548625461), (0.8611363115940526, 0.3478548451374538),
];
const GAUSS_LEGENDRE_8: [(f64, f64); 8] = [
    (-0.9602898564975363, 0.1012285362903763), (-0.7966664774136267, 0.2223810344533745),
    (-0.5255324099163290, 0.3137066458778873), (-0.1834346424956498, 0.3626837833783620),
    (0.1834346424956498, 0.3626837833783620), (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745), (0.9602898564975363, 0.1012285362903763),
];

/// Wrap a difference of two principal angles into (-pi, pi].
#[inline]
fn wrap(a: f64) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    if a > std::f64::consts::PI { a - two_pi } else if a <= -std::f64::consts::PI { a + two_pi } else { a }
}

/// DCT-II / DCT-III pair of length n, via an n-point complex FFT (Makhoul's algorithm).
///   forward:  X[k] = sum_j x[j] cos(pi k (2j+1) / 2n)
///   inverse:  recovers x from X.
struct Dct {
    n: usize,
    fwd: Arc<dyn Fft<f64>>,
    inv: Arc<dyn Fft<f64>>,
    twiddle: Vec<Complex<f64>>, // exp(-i pi k / 2n)
    buf: Vec<Complex<f64>>,
    scratch: Vec<Complex<f64>>,
}

impl Dct {
    fn new(n: usize) -> Dct {
        let mut planner = FftPlanner::new();
        let fwd = planner.plan_fft_forward(n);
        let inv = planner.plan_fft_inverse(n);
        let scratch_len = fwd.get_inplace_scratch_len().max(inv.get_inplace_scratch_len());
        let twiddle = (0..n).map(|k| Complex::from_polar(1.0, -std::f64::consts::PI * k as f64 / (2.0 * n as f64))).collect();
        Dct { n, fwd, inv, twiddle, buf: vec![Complex::new(0.0, 0.0); n], scratch: vec![Complex::new(0.0, 0.0); scratch_len] }
    }

    /// In place on x[offset + stride * j], j = 0..n
    fn forward(&mut self, x: &mut [f64], offset: usize, stride: usize) {
        let n = self.n;
        for j in 0..(n + 1) / 2 { self.buf[j] = Complex::new(x[offset + stride * 2 * j], 0.0); }
        for j in 0..n / 2 { self.buf[n - 1 - j] = Complex::new(x[offset + stride * (2 * j + 1)], 0.0); }
        self.fwd.process_with_scratch(&mut self.buf, &mut self.scratch);
        for k in 0..n { x[offset + stride * k] = (self.buf[k] * self.twiddle[k]).re; }
    }

    fn inverse(&mut self, x: &mut [f64], offset: usize, stride: usize) {
        let n = self.n;
        // V[k] = conj(twiddle[k]) (X[k] - i X[n-k]), X[n] = 0; v = ifft(V) / n
        for k in 0..n {
            let xk = x[offset + stride * k];
            let xnk = if k == 0 { 0.0 } else { x[offset + stride * (n - k)] };
            self.buf[k] = self.twiddle[k].conj() * Complex::new(xk, -xnk);
        }
        self.inv.process_with_scratch(&mut self.buf, &mut self.scratch);
        let s = 1.0 / n as f64;
        for j in 0..(n + 1) / 2 { x[offset + stride * 2 * j] = self.buf[j].re * s; }
        for j in 0..n / 2 { x[offset + stride * (2 * j + 1)] = self.buf[n - 1 - j].re * s; }
    }
}

/// Least-squares solver for the stream function on the corners of the physical region. Holds the
/// FFT plans, so keep one per geometry.
pub struct StreamFunctionSolver {
    m: usize,
    dct: Dct,
    eig: Vec<f64>, // eigenvalues of the weighted Neumann Laplacian, x-major
}

impl StreamFunctionSolver {
    pub fn new(geometry: &Geometry) -> StreamFunctionSolver {
        let nb = geometry.nboundary;
        let np = geometry.nx - 2 * nb;
        assert!(geometry.ny - 2 * nb == np, "stream function needs equal cell counts in x and y");
        let m = np + 1;
        let (wx, wy) = weights(geometry);
        // D^T D for the 1D Neumann chain of m nodes has DCT-II eigenvectors, eigenvalues 4 sin^2(pi k / 2m)
        let lam: Vec<f64> = (0..m).map(|k| 4.0 * (std::f64::consts::PI * k as f64 / (2.0 * m as f64)).sin().powi(2)).collect();
        let mut eig = vec![0.0; m * m];
        for a in 0..m { for b in 0..m { eig[a * m + b] = wx * lam[a] + wy * lam[b]; } }
        StreamFunctionSolver { m, dct: Dct::new(m), eig }
    }

    /// Solve (wx Dx^T Dx + wy Dy^T Dy) psi = rhs (x-major, m*m) in place, with mean(psi) = 0.
    fn solve(&mut self, rhs: &mut [f64]) {
        let m = self.m;
        for a in 0..m { self.dct.forward(rhs, a * m, 1); }
        for b in 0..m { self.dct.forward(rhs, b, m); }
        rhs[0] = 0.0;
        for k in 1..m * m { rhs[k] /= self.eig[k]; }
        for b in 0..m { self.dct.inverse(rhs, b, m); }
        for a in 0..m { self.dct.inverse(rhs, a * m, 1); }
    }
}

/// Least-squares weights of x- and y-segments, so that the sum approximates the integral of
/// |grad psi - g|^2 for non-square cells.
fn weights(geometry: &Geometry) -> (f64, f64) {
    let dx = geometry.delta_x();
    let dy = geometry.delta_y();
    (dy / dx, dx / dy)
}

/// The residual differences along x-segments (a, b) -> (a+1, b) and y-segments (a, b) -> (a, b+1)
/// of the corner grid (m*m each, x-major; the last x-segment of each row / y-segment of each
/// column is unused and zero).
pub(crate) fn residuals(geometry: &Geometry, elec_x: &Array2<f64>, elec_y: &Array2<f64>,
             flux_per_unit_charge: f64, xs: &[f64], ys: &[f64], qs: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let nb = geometry.nboundary;
    let m = geometry.nx - 2 * nb + 1;
    let dx = geometry.delta_x();
    let dy = geometry.delta_y();
    // psi_sing at the corners (principal values); differences along segments are wrapped
    let nq = qs.len();
    let mut theta = vec![0.0; nq * m * m];
    for k in 0..nq {
        for a in 0..m {
            for b in 0..m {
                theta[(k * m + a) * m + b] = (b as f64 * dy - ys[k]).atan2(a as f64 * dx - xs[k]);
            }
        }
    }
    let inv_two_pi = 1.0 / (2.0 * std::f64::consts::PI);
    // Close to a charge, the point charge's flux through a segment flips between +q/2 and -q/2 as the
    // charge crosses the segment's line, which would make psi_rest jump as a charge is dragged. So
    // the singular part subtracted here is regularised within a few cells, where the grid's charge
    // is smeared anyway: flux (q / 2 pi) * integral of (1 - exp(-r^2 / 2 s^2)) d theta, which is
    // continuous in the charge's position and equal to the point charge's beyond REG_FAR * s.
    let s_reg = REG_CELLS * dx.max(dy);
    let far2 = (REG_FAR * s_reg).powi(2);
    // near a charge: integrate the regularised d theta directly (16-point Gauss-Legendre)
    let regularised = |k: usize, x0: f64, y0: f64, x1: f64, y1: f64| -> f64 {
        let (ex, ey) = (x1 - x0, y1 - y0);
        let mut total = 0.0;
        for half in 0..2 {
            for (node, weight) in GAUSS_LEGENDRE_8.iter() {
                let t = 0.25 * (node + 1.0) + 0.5 * half as f64;
                let (rx, ry) = (x0 + t * ex - xs[k], y0 + t * ey - ys[k]);
                let r2 = rx * rx + ry * ry;
                // d theta / dt = (r x e) / r^2; regularised: (1 - exp(-r^2 / 2s^2)) / r^2 -> 1 / 2s^2
                let f = if r2 > 1e-12 * s_reg * s_reg { -(-r2 / (2.0 * s_reg * s_reg)).exp_m1() / r2 } else { 0.5 / (s_reg * s_reg) };
                total += 0.25 * weight * f * (rx * ey - ry * ex);
            }
        }
        total
    };
    // further away the segment cannot pass close to the charge: the exact (wrapped) d theta minus the
    // small, smooth correction integral of exp(-r^2 / 2s^2) d theta (4-point Gauss-Legendre)
    let correction = |k: usize, x0: f64, y0: f64, x1: f64, y1: f64| -> f64 {
        let (ex, ey) = (x1 - x0, y1 - y0);
        let mut total = 0.0;
        for (node, weight) in GAUSS_LEGENDRE_4.iter() {
            let t = 0.5 * (node + 1.0);
            let (rx, ry) = (x0 + t * ex - xs[k], y0 + t * ey - ys[k]);
            let r2 = rx * rx + ry * ry;
            total += 0.5 * weight * (-r2 / (2.0 * s_reg * s_reg)).exp() / r2 * (rx * ey - ry * ex);
        }
        total
    };
    let near2 = (REG_NEAR * s_reg).powi(2);
    let sing = |p0: usize, p1: usize| -> f64 {
        let (a0, b0, a1, b1) = (p0 / m, p0 % m, p1 / m, p1 % m);
        let (x0, y0, x1, y1) = (a0 as f64 * dx, b0 as f64 * dy, a1 as f64 * dx, b1 as f64 * dy);
        let (xm, ym) = (0.5 * (x0 + x1), 0.5 * (y0 + y1));
        let mut s = 0.0;
        for k in 0..nq {
            let d2 = (xm - xs[k]).powi(2) + (ym - ys[k]).powi(2);
            let dtheta = if d2 > far2 {
                wrap(theta[k * m * m + p1] - theta[k * m * m + p0])
            } else if d2 > near2 {
                wrap(theta[k * m * m + p1] - theta[k * m * m + p0]) - correction(k, x0, y0, x1, y1)
            } else {
                regularised(k, x0, y0, x1, y1)
            };
            s += qs[k] * dtheta;
        }
        s * inv_two_pi
    };
    let mut rx = vec![0.0; m * m];
    let mut ry = vec![0.0; m * m];
    let fx = dx / flux_per_unit_charge;
    let fy = dy / flux_per_unit_charge;
    for a in 0..m {
        for b in 0..m {
            let p = a * m + b;
            if a + 1 < m {
                // the face between cells (a, b-1) and (a, b) of the padded grid
                rx[p] = -elec_y[[a + nb, b + nb - 1]] * fx - sing(p, p + m);
            }
            if b + 1 < m {
                // the face between cells (a-1, b) and (a, b)
                ry[p] = elec_x[[a + nb - 1, b + nb]] * fy - sing(p, p + 1);
            }
        }
    }
    (rx, ry)
}

/// psi_rest on the (nphys+1)^2 corners of the physical region, x-major ([a][b] at a*(n+1)+b), with
/// mean zero.
pub fn stream_function_rest(solver: &mut StreamFunctionSolver, geometry: &Geometry, elec_x: &Array2<f64>, elec_y: &Array2<f64>,
                            flux_per_unit_charge: f64, xs: &[f64], ys: &[f64], qs: &[f64]) -> Vec<f64> {
    let m = solver.m;
    let (rx, ry) = residuals(geometry, elec_x, elec_y, flux_per_unit_charge, xs, ys, qs);
    let (wx, wy) = weights(geometry);
    // normal equations: rhs = wx Dx^T rx + wy Dy^T ry, (Dx^T r)[a] = r[a-1] - r[a] (r[-1] = r[m-1] = 0)
    let mut rhs = vec![0.0; m * m];
    for a in 0..m {
        for b in 0..m {
            let p = a * m + b;
            let mut v = 0.0;
            if a > 0 { v += wx * rx[p - m]; }
            if a + 1 < m { v -= wx * rx[p]; }
            if b > 0 { v += wy * ry[p - 1]; }
            if b + 1 < m { v -= wy * ry[p]; }
            rhs[p] = v;
        }
    }
    solver.solve(&mut rhs);
    rhs
}

// ---------------------------------------------------------------------------------------------
// Choice of the drawn levels (the gauge of psi)
//
// Field lines are the contours psi = k / LINES_PER_UNIT_CHARGE, so only psi modulo the line spacing
// matters, and adding a constant to psi shifts every line. Two choices:
//
// - canonical (static solvers, and the state a settled dynamic field relaxes to): psi_rest has mean
//   zero over the canvas (as solved). psi_rest is then odd under any reflection symmetry of the
//   configuration, so, as for the analytic models (decision D5), e.g. a +-q pair has a line along
//   its axis. (Matching psi to the free-space value around the edge instead does not work: the
//   periodic images make psi_rest vary there by most of a line spacing.)
// - causal (dynamic solver): away from currents, Maxwell's equations give
//       flux_per_unit_charge * d psi / dt = B_z + mean_jx * y - mean_jy * x
//   (the last terms from the uniform background current, see make_currents), exactly on the
//   staggered grid. So the offset is chosen each frame so that psi has changed by the time integral
//   of the right-hand side since the previous frame, at every corner away from the charges. Then psi
//   changes only where B does, and lines ahead of the wavefront do not move at all.
// ---------------------------------------------------------------------------------------------

/// Must match LINES_PER_UNIT_CHARGE in fieldmodels.js
pub const LINES_PER_UNIT_CHARGE: f64 = 4.0;

/// The drawn psi at the previous computation, for the causal choice.
pub struct GaugeHistory {
    /// psi_rest + offset + free-space psi of the point charges (principal values), at the corners
    pub psi: Vec<f64>,
    /// the charges' positions at that time
    pub positions: Vec<(f64, f64)>,
}

/// Free-space stream function of point charges (principal values) at the corners.
pub fn free_psi_at_corners(geometry: &Geometry, xs: &[f64], ys: &[f64], qs: &[f64]) -> Vec<f64> {
    let m = geometry.nx - 2 * geometry.nboundary + 1;
    let (dx, dy) = (geometry.delta_x(), geometry.delta_y());
    let mut out = vec![0.0; m * m];
    for a in 0..m {
        for b in 0..m {
            let (x, y) = (a as f64 * dx, b as f64 * dy);
            let mut s = 0.0;
            for k in 0..qs.len() { s += qs[k] * (y - ys[k]).atan2(x - xs[k]); }
            out[a * m + b] = s / (2.0 * std::f64::consts::PI);
        }
    }
    out
}

/// x - nearest integer
fn wrap_unit(x: f64) -> f64 { x - x.round() }

/// Circular mean (period 1) of values, robust against a few outliers: a second pass uses only the
/// values within 0.1 of the first estimate. None if there are too few values.
fn circular_mean(values: &[f64]) -> Option<f64> {
    if values.len() < 20 { return None; }
    let tau = 2.0 * std::f64::consts::PI;
    let mean_of = |vals: &mut dyn Iterator<Item = &f64>| -> Option<f64> {
        let (mut c, mut s, mut n) = (0.0, 0.0, 0);
        for v in vals { c += (tau * v).cos(); s += (tau * v).sin(); n += 1; }
        if n == 0 || (c * c + s * s).sqrt() < 1e-9 * n as f64 { None } else { Some(s.atan2(c) / tau) }
    };
    let first = mean_of(&mut values.iter())?;
    mean_of(&mut values.iter().filter(|v| wrap_unit(*v - first).abs() < 0.1)).or(Some(first))
}

/// Corners further than this from a charge (in cells, in each direction) are clear of its smeared
/// charge and current (4 sigma = 8 cells, plus a margin).
const EXCLUSION_CELLS: f64 = 10.0;

/// Is the corner (x, y) clear of the charges, including the paths along which they moved?
fn clear_of_charges(x: f64, y: f64, rx: f64, ry: f64, paths: &[((f64, f64), (f64, f64))]) -> bool {
    for &((x0, y0), (x1, y1)) in paths {
        // distance to the segment in units scaled by (rx, ry)
        let (px, py) = ((x - x0) / rx, (y - y0) / ry);
        let (ex, ey) = ((x1 - x0) / rx, (y1 - y0) / ry);
        let l2 = ex * ex + ey * ey;
        let t = if l2 > 0.0 { ((px * ex + py * ey) / l2).clamp(0.0, 1.0) } else { 0.0 };
        if (px - t * ex).powi(2) + (py - t * ey).powi(2) < 1.0 { return false; }
    }
    true
}

/// The causal offset (in units of charge) to add to psi_rest, given psi_raw = psi_rest + free psi
/// now, the history, the time integral of B_z (padded grid) and of the mean current since then, and
/// the charges' current positions. None if it can't be determined.
pub fn causal_offset(geometry: &Geometry, raw: &[f64], history: &GaugeHistory, b_integral: &Array2<f64>,
                     mean_current_integral: (f64, f64), flux_per_unit_charge: f64,
                     positions: &[(f64, f64)], deposit_positions: &[(f64, f64)]) -> Option<f64> {
    if history.positions.len() != positions.len() { return None; }
    let nb = geometry.nboundary;
    let m = geometry.nx - 2 * nb + 1;
    let (dx, dy) = (geometry.delta_x(), geometry.delta_y());
    let mut paths: Vec<_> = history.positions.iter().zip(positions.iter()).map(|(&p, &q)| (p, q)).collect();
    paths.extend(deposit_positions.iter().map(|&p| (p, p)));
    let (mx, my) = mean_current_integral;
    let mut vals = Vec::with_capacity(m * m);
    // the corners next to the edge are skipped: their faces border the absorbing layer
    for a in 2..m - 2 {
        for b in 2..m - 2 {
            let (x, y) = (a as f64 * dx, b as f64 * dy);
            if !clear_of_charges(x, y, EXCLUSION_CELLS * dx, EXCLUSION_CELLS * dy, &paths) { continue; }
            // psi at corner (a, b) changes with B_z at the Yee node of padded cell (a + nb - 1, b + nb - 1)
            let predicted = (b_integral[[a + nb - 1, b + nb - 1]] + mx * y - my * x) / flux_per_unit_charge;
            let p = a * m + b;
            vals.push(LINES_PER_UNIT_CHARGE * (history.psi[p] + predicted - raw[p]));
        }
    }
    circular_mean(&vals).map(|mu| mu / LINES_PER_UNIT_CHARGE)
}

/// The difference a - b of two offsets (units of charge), wrapped to within half a line spacing.
pub fn offset_difference(a: f64, b: f64) -> f64 {
    wrap_unit(LINES_PER_UNIT_CHARGE * (a - b)) / LINES_PER_UNIT_CHARGE
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fourier;

    fn geometry() -> Geometry { Geometry { x_max: 668.0, y_max: 335.0, nx: 111, ny: 111, nboundary: 13 } }

    fn fields() -> (Array2<f64>, Array2<f64>) {
        let mut ex = Array2::<f64>::zeros((111, 111));
        let mut ey = Array2::<f64>::zeros((111, 111));
        for i in 0..111 { for j in 0..111 {
            ex[[i, j]] = ((i * 7 + j * 3) % 11) as f64 * 10.0;
            ey[[i, j]] = ((i * 5 + j) % 13) as f64 * 10.0;
        } }
        (ex, ey)
    }

    /// Reference: the same least-squares problem solved by a periodic FFT of the mirrored grid.
    fn reference(g: &Geometry, rx: &[f64], ry: &[f64]) -> Vec<f64> {
        let m = g.nx - 2 * g.nboundary + 1;
        let n2 = 2 * m;
        let (wx, wy) = weights(g);
        let mirror = |e: usize| if e < m { e } else { n2 - 1 - e };
        let ext_x = |e: usize, f: usize| -> f64 {
            let fb = mirror(f);
            if e < m - 1 { rx[e * m + fb] } else if e == m - 1 || e == n2 - 1 { 0.0 } else { -rx[(n2 - 2 - e) * m + fb] }
        };
        let ext_y = |e: usize, f: usize| -> f64 {
            let ea = mirror(e);
            if f < m - 1 { ry[ea * m + f] } else if f == m - 1 || f == n2 - 1 { 0.0 } else { -ry[ea * m + n2 - 2 - f] }
        };
        let mut rhs = Array2::<Complex<f64>>::zeros((n2, n2));
        for e in 0..n2 { for f in 0..n2 {
            let em = (e + n2 - 1) % n2; let fm = (f + n2 - 1) % n2;
            rhs[[e, f]] = Complex::new(wx * (ext_x(em, f) - ext_x(e, f)) + wy * (ext_y(e, fm) - ext_y(e, f)), 0.0);
        } }
        let mut op = Array2::<Complex<f64>>::zeros((n2, n2));
        op[[0, 0]] = Complex::new(2.0 * wx + 2.0 * wy, 0.0);
        op[[1, 0]] = Complex::new(-wx, 0.0); op[[n2 - 1, 0]] = Complex::new(-wx, 0.0);
        op[[0, 1]] = Complex::new(-wy, 0.0); op[[0, n2 - 1]] = Complex::new(-wy, 0.0);
        fourier::array_fft(&mut op);
        fourier::array_fft(&mut rhs);
        rhs.zip_mut_with(&op, |r, o| { *r = if o.norm() > 1e-9 { *r / *o } else { Complex::new(0.0, 0.0) }; });
        fourier::array_ifft(&mut rhs);
        let mut out = vec![0.0; m * m];
        for a in 0..m { for b in 0..m { out[a * m + b] = rhs[[a, b]].re / (n2 * n2) as f64; } }
        let mean = out.iter().sum::<f64>() / (m * m) as f64;
        out.iter_mut().for_each(|v| *v -= mean);
        out
    }

    #[test]
    fn dct_roundtrip_and_definition() {
        for n in [1usize, 2, 5, 8, 86] {
            let mut d = Dct::new(n);
            let x: Vec<f64> = (0..n).map(|j| ((j * 37 + 11) % 17) as f64 - 8.0).collect();
            let mut y = x.clone();
            d.forward(&mut y, 0, 1);
            for k in 0..n {
                let want: f64 = (0..n).map(|j| x[j] * (std::f64::consts::PI * k as f64 * (2 * j + 1) as f64 / (2 * n) as f64).cos()).sum();
                assert!((y[k] - want).abs() < 1e-9 * (1.0 + want.abs()), "n={} k={} {} {}", n, k, y[k], want);
            }
            d.inverse(&mut y, 0, 1);
            for j in 0..n { assert!((y[j] - x[j]).abs() < 1e-9, "n={} j={}", n, j); }
        }
    }

    #[test]
    fn matches_mirrored_fft_reference() {
        let g = geometry();
        let (ex, ey) = fields();
        let xs = [100.0, 400.5]; let ys = [100.0, 200.25]; let qs = [2.0, -1.0];
        let mut solver = StreamFunctionSolver::new(&g);
        let psi = stream_function_rest(&mut solver, &g, &ex, &ey, 4000.0, &xs, &ys, &qs);
        let (rx, ry) = residuals(&g, &ex, &ey, 4000.0, &xs, &ys, &qs);
        let want = reference(&g, &rx, &ry);
        let scale = want.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        for (p, w) in psi.iter().zip(want.iter()) { assert!((p - w).abs() < 1e-9 * scale, "{} {}", p, w); }
    }

    #[test]
    fn exact_for_a_curl_free_divergence_free_flux() {
        // a uniform field has a stream function: psi_rest must reproduce it exactly
        let g = geometry();
        let ex = Array2::<f64>::from_elem((111, 111), 30.0);
        let ey = Array2::<f64>::from_elem((111, 111), -20.0);
        let mut solver = StreamFunctionSolver::new(&g);
        let psi = stream_function_rest(&mut solver, &g, &ex, &ey, 4000.0, &[], &[], &[]);
        let m = 86;
        let (dx, dy) = (g.delta_x(), g.delta_y());
        // psi = (Ex y - Ey x) / flux_per_unit_charge + const
        let f = |a: usize, b: usize| (30.0 * b as f64 * dy + 20.0 * a as f64 * dx) / 4000.0;
        let mean = (0..m * m).map(|p| f(p / m, p % m)).sum::<f64>() / (m * m) as f64;
        for p in 0..m * m { assert!((psi[p] - (f(p / m, p % m) - mean)).abs() < 1e-9, "{}", p); }
    }
}
