use crate::Geometry;

pub struct PmlIterator {
    i: usize,
    j: usize,
    nx: usize,
    ny: usize,
    layer_width: usize,
}

impl Iterator for PmlIterator {
    type Item = (usize, usize, f64, f64); // (i, j, sigma_x, sigma_y)

    fn next(&mut self) -> Option<Self::Item> {
        self.i+=1;
        if self.i == self.nx {
            self.i = 0;
            self.j += 1;
            if self.j == self.ny {
                return None;
            }
        }
        let sigma_x = if self.i < self.layer_width || self.i >= self.nx - self.layer_width { 
            let distance_to_boundary = self.i.min(self.nx - self.i - 1) as f64;
            1.0 - distance_to_boundary / self.layer_width as f64 
        } else { 0.0 };
        let sigma_y = if self.j < self.layer_width || self.j >= self.ny - self.layer_width { 
            let distance_to_boundary = self.j.min(self.ny - self.j - 1) as f64;
            1.0 - distance_to_boundary / self.layer_width as f64 
        } else { 0.0 };
        Some((self.i, self.j, sigma_x, sigma_y))
    }
}

pub fn pml_iterator_from_geometry(geometry: &Geometry) -> PmlIterator {
    PmlIterator { i: 0, j: 0, nx: geometry.nx, ny: geometry.ny, layer_width: geometry.nboundary }
}