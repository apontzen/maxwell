#[derive(Clone)]
pub struct Geometry {
    /// maximum x extent of the physical region. Physical region runs from 0 to x_max; boundary cells extend further.
    pub x_max: f64, 
    /// maximum y extent of the physical region. Physical region runs from 0 to y_max; boundary cells extend further.
    pub y_max: f64, 
    /// number of grid cells in x direction, including boundary cells (=2*nboundary)
    pub nx: usize, 
    /// number of grid cells in y direction, including boundary cells (=2*nboundary)
    pub ny: usize, 
    /// number of boundary cells on each side
    pub nboundary: usize, 
}

impl Geometry {
    pub fn delta_x(&self) -> f64 {
        self.x_max / (self.nx - 2 * self.nboundary) as f64
    }

    pub fn delta_y(&self) -> f64 {
        self.y_max / (self.ny - 2 * self.nboundary) as f64
    }

    pub fn position_to_cell(&self, x: f64, y: f64) -> Option<(usize, usize)> {
        let (i, j) = self.position_to_cell_unclamped(x, y);
        if i < 0 || i >= self.nx as isize || j < 0 || j >= self.ny as isize {
            None
        } else {
            Some((i as usize, j as usize))
        }
    }

    pub fn position_to_surrounding_cells(&self, x: f64, y: f64) -> Vec<(usize, usize)> {
        let (i, j) = self.position_to_cell_unclamped(x, y);
        let mut result = vec![];
        for i_offset in -1..2 {
            for j_offset in -1..2 {
                let i = i + i_offset;
                let j = j + j_offset;
                if i >= 0 && i < self.nx as isize && j >= 0 && j < self.ny as isize {
                    result.push((i as usize, j as usize));
                }
            }
        }
        result
    }

    pub fn position_to_cell_unclamped(&self, x: f64, y: f64) -> (isize, isize) {
        let i = (x / self.x_max * (self.nx - 2*self.nboundary) as f64) as isize + self.nboundary as isize;
        let j = (y / self.y_max * (self.ny - 2*self.nboundary) as f64) as isize + self.nboundary as isize;
        (i, j)
    }

    pub fn cell_to_centroid(&self, i: usize, j: usize) -> (f64, f64) {
        let x = (i as f64 - self.nboundary as f64 + 0.5) * self.delta_x();
        let y = (j as f64 - self.nboundary as f64 + 0.5) * self.delta_y();
        (x, y)
    }

    pub fn cell_to_corners(&self, i: usize, j: usize) -> Vec<(f64, f64)> {
        let x0 = (i as f64 - self.nboundary as f64) * self.delta_x();
        let x1 = (i as f64 - self.nboundary as f64 + 1.0) * self.delta_x();
        let y0 = (j as f64 - self.nboundary as f64) * self.delta_y();
        let y1 = (j as f64 - self.nboundary as f64 + 1.0) * self.delta_y();
        
        vec![(x0, y0), (x1, y0), (x0, y1), (x1, y1)]
    }

    pub fn x_extent_including_boundary(&self) -> f64 {
        self.x_max * self.nx as f64 / (self.nx - 2*self.nboundary) as f64
    }

    pub fn y_extent_including_boundary(&self) -> f64 {
        self.y_max * self.ny as f64 / (self.ny - 2*self.nboundary) as f64
    }

    pub fn in_padding_region(&self, x: f64, y: f64) -> bool {
        let min_x = -(self.nboundary as f64) * self.delta_x();
        let max_x = (self.nx as f64 - self.nboundary as f64) * self.delta_x();
        let min_y = -(self.nboundary as f64) * self.delta_y();
        let max_y = (self.ny as f64 - self.nboundary as f64) * self.delta_y();
        x>=min_x && x<=max_x && y>=min_y && y<=max_y
    }

}
