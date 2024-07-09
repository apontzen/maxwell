use wasm_bindgen::prelude::*;
#[allow(unused)]
use serde_wasm_bindgen::{from_value, to_value};
use ndarray::Array2;
use crate::{Pair, FieldConfiguration};
use web_sys::console;

type PotentialCalculator = fn(&FieldConfiguration, f64, f64) -> f64;
type PotentialGradientCalculator = fn(&FieldConfiguration, f64, f64) -> Pair;

struct ContouringCollection<'a> {
    potential_calculator: PotentialCalculator,
    potential_gradient_calculator: PotentialGradientCalculator,
    configuration: &'a FieldConfiguration
}

/// Generate a contour of the potential field starting at the specified point
fn generate_potential_contours(description: &ContouringCollection, x0: f64, y0: f64, 
    level: f64, reverse_direction: bool) -> Vec<(f64, f64)> {
    
    let step_size: f64 = if reverse_direction { -0.5 } else { 0.5 };

    let pot_grad = description.potential_gradient_calculator;
    let field_configuration = description.configuration;


    const RETURN_EVERY: i32 = 10;
    const FINISH_TOLERANCE: f64 = 0.5;
    const MINIMUM_NUM_STEPS: i32 = 20;
    const MAXIMUM_NUM_STEPS: i32 = 10000; // just to prevent hangs in weird situations

    let mut x = x0;
    let mut y = y0;
    let mut contour: Vec<(f64, f64)> = vec![(x, y)];
    let mut step = 0;


    loop {
        // predictor step...
        let Pair {u, v} = pot_grad(field_configuration, x, y);
        let r = (u * u + v * v).sqrt();
        
        x -= step_size * v / r;
        y += step_size * u / r;

        // corrector step...
        let Pair {u:u1, v: v1} = pot_grad(field_configuration, x, y);
        let r1 = (u1 * u1 + v1 * v1).sqrt();

        x -= step_size * 0.5 * (v1 - v) / r1;
        y += step_size * 0.5 * (u1 - u) / r1;

        step += 1;

        if step % RETURN_EVERY == 0 {
            contour.push((x, y));
        }
        if !field_configuration.geometry.in_padding_region(x, y) {
            if !reverse_direction {
                // If the contour doesn't close in on itself, we need to supplement it with the contour
                // running in the opposite direction from the start point
                let contour_reversed = generate_potential_contours(description, x0, y0, level, true);

                // now add contours_reversed onto the beginning of contour, in reverse order
                contour_reversed.iter().for_each(|(x, y)| contour.insert(0, (*x, *y))); 
            }

            break;
        }
        if step >= MAXIMUM_NUM_STEPS {
            break;
        }

        if ((x-x0).powi(2) + (y-y0).powi(2)) < FINISH_TOLERANCE.powi(2) && step > MINIMUM_NUM_STEPS {
            break;
        }
    }
    contour.push((x, y));
    contour
}


/// Find a nearby point to (x0, y0) at which the potential has the specified level
fn find_crossing_point(description: &ContouringCollection, level: f64, x0: f64, y0: f64) -> (f64, f64) {
    let mut x = x0;
    let mut y = y0;

    let field_configuration = description.configuration;
    let pot_calc = description.potential_calculator;
    let pot_grad = description.potential_gradient_calculator;


    let mut potential = pot_calc(field_configuration, x, y);
    let mut step = 0;

    const MAXIMUM_NUM_STEPS: i32 = 200;
    const TOLERANCE: f64 = 1e-4;
    
    let orig_potential = potential;
    

    if (potential-level).abs() < TOLERANCE {
        return (x0, y0);
    }
    let potential_offset_sign = (potential-level).signum();

    

    while (potential-level).signum() == potential_offset_sign {
        let Pair {u, v} = pot_grad(field_configuration, x, y);
        let r = (u * u + v * v).sqrt();
        let (u_normed, v_normed) = (u/r.powi(2), v/r.powi(2));

        if r<1e-6 || step>MAXIMUM_NUM_STEPS {
            console::log_1(&format!("Searching for crossing point from ({}, {}) at level {}; starting {}", x0, y0, level, orig_potential).into());
            console::log_1(&format!("Aborted after {} steps @ {} {} : {}", step, x, y, potential).into());

            return (x0, y0);
        }
        // try to overshoot the crossing point, then we'll bisect to find it exactly
        x -= 1.5*(potential-level)*u_normed;
        y -= 1.5*(potential-level)*v_normed;
        potential = pot_calc(field_configuration, x, y);
        step += 1;
       
    } 
    let (mut x0, mut y0) = (x0, y0);
    let (mut x1, mut y1) = (x, y);    

    // now do bisection search between (x0, y0) and (x1, y1) to find crossing point
    step = 0;
   
    while (potential-level).abs() > TOLERANCE {
        x = (x0 + x1) / 2.0;
        y = (y0 + y1) / 2.0;
        potential = pot_calc(field_configuration, x, y);

        if (potential-level).signum() == potential_offset_sign {
            x0 = x;
            y0 = y;
        } else {
            x1 = x;
            y1 = y;
        }
        
        step += 1;
        if step > MAXIMUM_NUM_STEPS {
            let final_pot = pot_calc(field_configuration, x, y);
            console::log_1(&format!("Aborted bisection ({} {}): {} point after {} steps", x, y, final_pot, step).into());
            break;
        }
    }
    
    

    (x,y)
}

#[wasm_bindgen]
pub fn generate_potential_contours_at_level(field_configuration: &FieldConfiguration, level: f64) -> JsValue {
    /// Generate a contour at a specified level of the potential field
        
    const MAX_CONTOURS: usize = 5;

    let description  = ContouringCollection { potential_calculator: crate::compute_potential_electrostatic_direct,
                                                            potential_gradient_calculator: crate::compute_field_electrostatic_direct,
                                                            configuration: field_configuration};
    
    // First, scan over the grid to find all cells that cross the level
    let mut crossing_flags = Array2::<bool>::from_shape_fn((field_configuration.geometry.nx, field_configuration.geometry.ny), |(i, j)| {
        let change_in_cell = field_configuration.geometry.cell_to_corners(i, j).
            iter().
            any(|(x, y)| (description.potential_calculator)(field_configuration, *x, *y) > level) !=
            field_configuration.geometry.cell_to_corners(i, j).
            iter().
            all(|(x, y)| (description.potential_calculator)(field_configuration, *x, *y) > level);

        change_in_cell
    });

    let mut contours : Vec<Vec<(f64, f64)>> = vec![];

    loop {
        if contours.len() >= MAX_CONTOURS {
            console::log_1(&format!("Ran out of contours for {}",level).into());
            break;
        }


        // Find the first cell that crosses the level
        let (i, j) = match (0..field_configuration.geometry.nx).flat_map(|i| (0..field_configuration.geometry.ny).map(move |j| (i, j))).find(|(i, j)| crossing_flags[[*i, *j]]) {
            Some((i, j)) => (i, j),
            None => break,
        };

        crossing_flags[[i, j]] = false;

        // Now follow the contour from this cell
        let (mut x, mut y) = field_configuration.geometry.cell_to_centroid(i, j);
        (x, y) = find_crossing_point(&description, level, x, y);
        
        
        let contour = generate_potential_contours(&description, x, y, level, false);

        // unflag cells that have been visited by this contour
        contour.iter().for_each(|(x, y)| {
            field_configuration.geometry.position_to_surrounding_cells(*x, *y).iter().for_each(
                |(it, jt)| { crossing_flags[[*it,*jt]] = false; }
            );
            /*match field_configuration.geometry.position_to_cell(*x, *y) {
                Some((it, jt)) => crossing_flags[[it, jt]] = false,
                None => (),
            }*/
        });


        contours.push(contour);


    }

    // web_sys::console::log_1(&format!("Found {} contours", contours.len()).into());

    to_value(&contours).unwrap()


}