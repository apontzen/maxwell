use wasm_bindgen::prelude::*;
#[allow(unused)]
use serde_wasm_bindgen::{from_value, to_value};
use ndarray::Array2;
use crate::{Pair, FieldConfiguration, Charge};
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

        if ((description.potential_calculator)(field_configuration, x, y) - level).abs() > 0.1 {
            (x,y) = find_crossing_point(description, level, x, y);
        }

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
        // {u,v}_normed points down the potential gradient, and we want to move down if potential>level, up if potential<level
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

fn squared_distance_to_closest_charge(field_configuration: &FieldConfiguration, x: f64, y: f64) -> f64 {
    let closest_charge = field_configuration.closest_charge(x,y).unwrap(); // this definitely exists
    (x-closest_charge.x).powi(2) + (y-closest_charge.y).powi(2)
}

fn generate_contours_at_levels(description: &ContouringCollection, levels: Vec<f64>) -> Vec<Vec<(f64, f64)>> {
    const MAX_CONTOURS: usize = 50;
    const NO_CHANGE: f64 = f64::INFINITY;

    let field_configuration = description.configuration;

    
    // First, scan over the grid to find all cells that cross the level
    let mut crossing_level = Array2::<f64>::from_shape_fn((field_configuration.geometry.nx, field_configuration.geometry.ny), |(i, j)| {
        let potential_in_corners = field_configuration.geometry.cell_to_corners(i, j).
            iter().
            map(|(x, y)| (description.potential_calculator)(field_configuration, *x, *y)).
            collect::<Vec<f64>>();

        for level in &levels {
            let change_in_cell = potential_in_corners.iter().
                any(|&potential| potential > *level) !=
                potential_in_corners.iter().
                all(|&potential| potential > *level);
                
            if change_in_cell {
                return *level;
            }
        }
        NO_CHANGE
    });

    let mut contours : Vec<Vec<(f64, f64)>> = vec![];

    loop {
        if contours.len() >= MAX_CONTOURS {
            console::log_1(&format!("Ran out of contours").into());
            break;
        }


        // Find the first cell that crosses the level
        let (i, j) = match (0..field_configuration.geometry.nx).
             flat_map(|i| (0..field_configuration.geometry.ny).
             map(move |j| (i, j))).
             find(|(i, j)| crossing_level[[*i, *j]]!=NO_CHANGE) 
             {
                Some((i, j)) => (i, j),
                None => break,
             };

        let level = crossing_level[[i, j]];

        crossing_level[[i, j]] = NO_CHANGE;

        // Now follow the contour from this cell
        let (mut x, mut y) = field_configuration.geometry.cell_to_centroid(i, j);

        if squared_distance_to_closest_charge(&field_configuration, x, y) < 100. {
            continue;
        }


        
        (x, y) = find_crossing_point(&description, level, x, y);

        
        
        let contour = generate_potential_contours(&description, x, y, level, false);

        let mut max_squared_distance = 0.0;

        // unflag cells that have been visited by this contour
        contour.iter().for_each(|(x, y)| {
            let squared_distance = squared_distance_to_closest_charge(&field_configuration, *x, *y);
            if squared_distance > max_squared_distance {
                max_squared_distance = squared_distance;
            }

            field_configuration.geometry.position_to_surrounding_cells(*x, *y).iter().for_each(
                |(it, jt)| { 
                    if crossing_level[[*it,*jt]]==level { 
                        crossing_level[[*it,*jt]] = NO_CHANGE; 
                    }
                }
            );
        });


        if max_squared_distance > 150.0 {
            // only add the contour if it's not too close to a charge
            contours.push(contour);
        }

    }
    contours 

}


/// Generate a contour at a specified level of the electrostatic potential field
#[wasm_bindgen]
pub fn generate_potential_contours_at_levels(field_configuration: &FieldConfiguration, levels: Vec<f64>) -> JsValue {
    let description  = ContouringCollection { potential_calculator: crate::compute_potential_electrostatic_direct,
        potential_gradient_calculator: crate::compute_field_electrostatic_direct,
        configuration: field_configuration};

    let contours = generate_contours_at_levels(&description, levels);

    to_value(&contours).unwrap()
}

fn line_crosses_symmetry(field_configuration: &FieldConfiguration, x0: f64, y0: f64, x1: f64, y1: f64) -> Option<Pair> {

    let crosses_symmetry = |charge: &crate::Charge, other_charge: &crate::Charge| -> Option<Pair> {
        let f = |x,y| {(x-charge.x)*(other_charge.y-charge.y) - (y-charge.y)*(other_charge.x-charge.x)};
        let f0 = f(x0, y0);
        let f1 = f(x1, y1);
        if f0.signum() != f1.signum() {
            let t = f0 / (f0 - f1);
            let x = x0 + t * (x1 - x0);
            let y = y0 + t * (y1 - y0);
            Some(Pair {u: x, v: y})
        } else {
            None
        }
    };

    if field_configuration.charges.len() == 1 {
        // pretend there is another charge displaced to the right of the one charge
        let charge = &field_configuration.charges[0];
        let other_charge = crate::Charge {x: charge.x + 1.0, y: charge.y, charge: charge.charge};
        return crosses_symmetry(charge, &other_charge);
    }

    let score_pair = |charge: &crate::Charge, other_charge: &crate::Charge| -> f64 {
        let r_squared = ((charge.x - other_charge.x).powi(2) + (charge.y - other_charge.y).powi(2));
        (1.3-charge.charge*other_charge.charge)/r_squared
    };

    let mut charge_already_paired: Vec<bool> = vec![false; field_configuration.charges.len()];

    for i in 0..field_configuration.charges.len() {
        if charge_already_paired[i]
        {
            continue;
        }

        let charge = &field_configuration.charges[i];
        // find the highest-scoring other_charge:
        let mut best_score = 0.0;
        let mut best_other_charge = None;
        let mut best_other_charge_index = 0;
        for j in 0..field_configuration.charges.len() {
            let other_charge = &field_configuration.charges[j];
            if charge == other_charge || charge_already_paired[j] {
                continue;
            }
            let score = score_pair(charge, other_charge);
            if score > best_score {
                best_score = score;
                best_other_charge = Some(other_charge);
                best_other_charge_index = j;
            }
        }

        if best_other_charge.is_none() {
            continue;
        }

        charge_already_paired[i] = true;
        charge_already_paired[best_other_charge_index] = true;
        
        match crosses_symmetry(charge, best_other_charge.unwrap()) {
            Some(pair) => return Some(pair),
            None => (),
        }
            
        
    }

    None
}

#[wasm_bindgen]
pub fn generate_potential_contours_and_arrow_positions_at_levels(field_configuration: &FieldConfiguration, levels: Vec<f64>) -> JsValue {
    let description  = ContouringCollection { potential_calculator: crate::compute_potential_electrostatic_direct,
        potential_gradient_calculator: crate::compute_field_electrostatic_direct,
        configuration: field_configuration};

    let contours = generate_contours_at_levels(&description, levels);
    let mut arrows: Vec<(f64, f64)> = vec![];

    let mut steps_until_another_arrow_allowed: usize = 0;
    
    for contour in &contours {
        for ((x0, y0), (x1, y1)) in contour.iter().zip(contour.iter().skip(1)) {
            if steps_until_another_arrow_allowed > 0 {
                steps_until_another_arrow_allowed -= 1;
                continue;
            }
            match line_crosses_symmetry(field_configuration, *x0, *y0, *x1, *y1) {
                Some(Pair {u, v}) => { 
                    steps_until_another_arrow_allowed = 10;
                    arrows.push((u, v))
                },
                None => (),
            }
        }
    }

    to_value(&(contours, arrows)).unwrap()
}
