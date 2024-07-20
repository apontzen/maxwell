import { compute_field_electrostatic_direct, generate_potential_contours_at_levels, generate_potential_contours_and_arrow_positions_at_levels, compute_field_magnetostatic_direct, compute_field_electrostatic_direct_to_buffer } from './maxwell/out/maxwell.js';
import { getChargeFromPoint, chargeSize } from './draw.js';

const DEBUG_MESSAGES = false; // Warning: can generate a LOT of console output!

function debug_log(...message) {
    if(DEBUG_MESSAGES) {
        console.log(...message);
    }
}
class MultiStreamDepartures {
    constructor(num_departures, starting_angle = 0) {
        this.num_departures = num_departures;
        this.next_departure_index = 0;
        this.starting_angle = starting_angle;
        this.arrival_range = {};
        this.arrivals = {}
        this.departure_ranges = null;

    }

    get_new_departure_no_arrivals(departure_index_interleaved) {
        return this.starting_angle + (departure_index_interleaved+0.5) * 2*Math.PI / this.num_departures;
    }

    calculate_departure_ranges() {
        const TOL = 1e-5; // a 'small' tolerance angle, used to avoid numerical issues

        // now the space available for actually launching streamlines is the space left after
        // excluding all arrival ranges. To complicate matters, the arrival_ranges may overlap, and
        // may also wrap around the 2 pi interval. So we need to find all the spaces left after 
        // excluding the arrival ranges, with these two complications in mind.
        let available_ranges = [{start: 0, end: 2*Math.PI}];
        let exclusion_ranges = []

        // step 0: to try and maintain rotational invariance, orient the available range so that it
        // starts at the first known arrival
        if (Object.entries(this.arrival_range).length === 1) {
            let first_arrival = Math.min(...Object.values(this.arrival_range).map(([start, end]) => start));
            available_ranges = [{start:0, end: first_arrival-TOL}, {start: first_arrival+TOL, end: first_arrival + 2*Math.PI}];
        }

        // step 1: if any exclusion range wraps around, turn it into two exclusion ranges
        for (let [, [exclusion_start, exclusion_end]] of Object.entries(this.arrival_range)) {
            if (exclusion_start<=exclusion_end) {
                exclusion_ranges.push([exclusion_start, exclusion_end]);
            } else {
                exclusion_ranges.push([0, exclusion_end]);
                exclusion_ranges.push([exclusion_start, 2*Math.PI]);
            }
        }

        // step 2: iteratively break launch ranges into smaller ranges by excluding the exclusion ranges
        for (let [exclusion_start, exclusion_end] of exclusion_ranges) {
            if (exclusion_start == exclusion_end) {
                continue;
            }
            let new_ranges = [];
            for (let range of available_ranges) { 
                if (range.start<=exclusion_start && range.end>=exclusion_start) {
                    new_ranges.push({start: range.start, end: exclusion_start});
                }
                if (range.end>=exclusion_end && range.start<=exclusion_end) {
                    new_ranges.push({start: exclusion_end, end: range.end});
                }
                if ((exclusion_start<=range.start && exclusion_end<=range.start) 
                    || (exclusion_start>=range.end && exclusion_end>=range.end)) {
                    new_ranges.push(range);
                }
            }
            available_ranges = new_ranges;
        }

        debug_log("departure ranges (end step 2):", available_ranges);
        
        // step 3: join up any range that starts at 0 with any range that ends at 2 pi
        if(available_ranges.length>1) {
            let new_ranges = [];
            
            for (let range of available_ranges) {
                if (range.start <= TOL) {
                    for (let range2 of available_ranges) {
                        if (range2.end >= 2*Math.PI - TOL) {
                            new_ranges.push({start: range2.start, end: range.end + 2*Math.PI});
                        }
                    }
                } else if (range.end < 2*Math.PI - TOL) {
                    new_ranges.push(range);
                }
            }

            available_ranges = new_ranges;
        }
        this.departure_ranges = available_ranges;


        let total_available_range = 0;
        for (let range of available_ranges) {
            total_available_range += range.end - range.start;
        }

        let num_departure_per_range = available_ranges.map(
            range => Math.round(this.num_departures*(range.end - range.start)/total_available_range)
        );

        // if sum(num_departure_per_range) is less than num_departures, add one to the first
        // range 
        let num_departure_per_range_sum = num_departure_per_range.reduce((a, b) => a + b, 0);
        if (num_departure_per_range_sum < this.num_departures) {
            num_departure_per_range[0]+=this.num_departures - num_departure_per_range_sum;
        }

        this.num_departure_per_range = num_departure_per_range;

        debug_log("arrivals:", this.arrivals);
        debug_log("exclusion ranges:", exclusion_ranges);
        debug_log("departure ranges:", available_ranges);
        debug_log("num_departure_per_range:", num_departure_per_range);

    }

    get_new_departure_with_arrivals(departure_index_interleaved) {
        let logical_angle = (departure_index_interleaved+0.5)/this.num_departures;
        if (this.departure_ranges === null) {
            this.calculate_departure_ranges();
        }

        let available_ranges = this.departure_ranges;
        let num_departure_per_range = this.num_departure_per_range;
        

        let index_within_range = departure_index_interleaved;
        let range_index = 0;
        while (index_within_range >= num_departure_per_range[range_index]) {
            index_within_range -= num_departure_per_range[range_index];
            range_index++;
        }
        
        let range_min = available_ranges[range_index].start;
        let range_max = available_ranges[range_index].end;
        let angle = range_min + (range_max - range_min) * (1 + index_within_range) / (1 + num_departure_per_range[range_index]);

        return angle;
 
    }

    get_new_departure() {
        const departure_index = this.next_departure_index;
        if(departure_index>=this.num_departures)
            return null;
        let angle; 

        // in case we find there is a problem with a streamline at one end of our range, we
        // interleave the departures, so that while departure_index runs from 0 to n-1 inclusive,
        // departure_index_interleaved jumps from 0 to n-1 to 1 to n-2 to 2 to n-3, etc.
        // Then, if a problem is found at one end of the range, the correction is made and the
        // streamlines don't bunch up around that end.
        const departure_index_interleaved = (departure_index%2==0)?(departure_index/2):(this.num_departures - 1 - (departure_index-1)/2);

        const restricted_angle = Object.keys(this.arrival_range).length > 0;

        if (restricted_angle) {
            angle = this.get_new_departure_with_arrivals(departure_index_interleaved);
        } else {
            angle = this.get_new_departure_no_arrivals(departure_index_interleaved);
        }

        this.next_departure_index++;

        return angle;
    }

    arrivals_to_range(arrivals, source_charge, arrival_charge) {

        const angle_from_source = Math.atan2(source_charge.y - arrival_charge.y, source_charge.x - arrival_charge.x);

        let min_angle_relative_to_source = 2*Math.PI;
        let max_angle_relative_to_source = -2*Math.PI;

        for (let arrival of arrivals) {
            let angle_relative_to_source = arrival - angle_from_source;
            if (angle_relative_to_source > Math.PI) {
                angle_relative_to_source -= 2*Math.PI;
            }
            if (angle_relative_to_source < -Math.PI) {
                angle_relative_to_source += 2*Math.PI;
            }
            if (angle_relative_to_source < min_angle_relative_to_source) {
                min_angle_relative_to_source = angle_relative_to_source;
            }
            if (angle_relative_to_source > max_angle_relative_to_source) {
                max_angle_relative_to_source = angle_relative_to_source;
            }
        }

        let min_angle = angle_from_source + min_angle_relative_to_source;
        let max_angle = angle_from_source + max_angle_relative_to_source;


        if(min_angle<0)
            min_angle+=2*Math.PI;

        if(max_angle<0)
            max_angle+=2*Math.PI;

        if(min_angle>2*Math.PI)
            min_angle-=2*Math.PI;

        if(max_angle>2*Math.PI)
            max_angle-=2*Math.PI;

        debug_log("arrivals:", arrivals, "min_angle:", min_angle, "max_angle:", max_angle, "angle_from_source", angle_from_source);


        return [min_angle, max_angle];
    }

    register_arrival(angle, source_charge, arrival_charge) {

        const source_id = source_charge.id;

        angle = (angle + 2*Math.PI)%(2*Math.PI);

        this.num_departures--;
        if (!(source_id in this.arrivals)) {
            this.arrivals[source_id] = [];
        }
        this.arrivals[source_id].push(angle);

        this.arrival_range[source_id] = this.arrivals_to_range(this.arrivals[source_id], source_charge, arrival_charge);
        this.departure_ranges = null; // important: invalidate departure ranges
        
    }
}

function initialiseFieldlineDrawingInfo(charges) {
    // Decide which direction the first field line should be launched from a charge
    // We want this to be rotationally invariant, translationally invariant, and invariant to the
    // ordering of the charges in the list. We also want it to reflect any symmetry in the charge
    // distribution. We can achieve this by launching the first field line towards the centre of
    // mass of the charge. We weight each charge by its absolute charge magnitude, so that
    // large charges have more influence on the direction of the field line than small charges.
    const numCharges = charges.length;
    if (numCharges === 1) {
        charges[0].angle = 0;
        charges[0].score = 1;
        charges[0].processed = false;
        charges[0].departures = new MultiStreamDepartures(Math.abs(charges[0].charge) * 4, 0);
    } else if (numCharges > 1) {
        let com_x = 0;
        let com_y = 0;
        let total_mass = 0;
        for (let i = 0; i < numCharges; i++) {
            const charge = charges[i];
            com_x += Math.abs(charge.charge) * charge.x;
            com_y += Math.abs(charge.charge) * charge.y;
            total_mass += Math.abs(charge.charge);
        }

        // Unlikely exceptional case: if all charges are zero (?!), just pick the first charge
        if (total_mass === 0) {
            com_x = charges[0].x;
            com_y = charges[0].y;
        } else {
            com_x /= total_mass;
            com_y /= total_mass;
        }

        for (let i = 0; i < numCharges; i++) {
            const currentCharge = charges[i];
            const dx = com_x - currentCharge.x;
            const dy = com_y - currentCharge.y;
            if (dx === 0 && dy === 0) {
                currentCharge.angle = 0;
            } else {
                currentCharge.angle = Math.atan2(dy, dx);
            }
            currentCharge.score = Math.sqrt(dx * dx + dy * dy);
            currentCharge.processed = false;
            currentCharge.departures = new MultiStreamDepartures(Math.abs(currentCharge.charge) * 4, currentCharge.angle);
        }
        
    }

}

function clearFieldlineDrawingInfo(charges) {
    for (let charge of charges) {
        delete charge.angle;
        delete charge.score;
        delete charge.processed;
        delete charge.departures;
    }
}

function getNextChargeToProcess(charges) {
    let highest_score = -1;
    let highest_score_charge = null;
    for (let charge of charges) {
        const score = charge.score; //  + Object.keys(charge.departures.arrivals).length * 1000;
        if (!charge.processed && score > highest_score) {
            highest_score = score;
            highest_score_charge = charge;
        }
    }
    if(highest_score_charge!==null)
        highest_score_charge.processed = true;

    return highest_score_charge;
}


export function drawElectrostaticFieldLines(charges, field, ctx, rect, chargeSize) {

    initialiseFieldlineDrawingInfo(charges);

    // at each step, a score is calculated for each charge, and the highest-scoring charge is processed
    // next. To start with, the score is the distance of the charge from the "centre of the field" (defined in
    // fieldlineStartingAngles).
    // 
    // However, if a field line arrives at a charge, processing that charge becomes more urgent, so 1000 is 
    // added to the score of that charge. This is to try and distribute field lines evenly around that charge
    // rather than risk launching more fieldlines at it and getting very uneven coverage.

    let charge;

    while(charge = getNextChargeToProcess(charges)) {
    
        const x = charge.x;
        const y = charge.y;
        const departures = charge.departures;
        
        debug_log("PROCESS CHARGE:",charge);

        let num_failed_launches = 0;
        const max_failed_launches = 20;

        while(true) {
            
            const stream_angle = departures.get_new_departure();
            if (stream_angle === null) break;

            debug_log("Departure number", departures.next_departure_index, "of", departures.num_departures, "angle:", stream_angle);

            let stream_x = x + chargeSize * Math.cos(stream_angle);
            let stream_y = y + chargeSize * Math.sin(stream_angle);
            let length_covered = 0;

            
            

            
            const step_size = 1.0;
            const step = charge.charge>0?step_size:-step_size;

            const x_steps = [stream_x];
            const y_steps = [stream_y];

            let n_steps = 0;

            let u_last = 0;
            let v_last = 0;

            let buffer = new Float64Array(2);

            while((length_covered<20 || getChargeFromPoint(charges, stream_x, stream_y) === null)
                && (stream_x>-rect.width && stream_y>-rect.height && stream_x<2*rect.width && stream_y<2*rect.height)
                && n_steps<5000) {
                n_steps++;

                compute_field_electrostatic_direct_to_buffer(field, stream_x, stream_y, buffer);
                let u = buffer[0];
                let v = buffer[1];
                let norm = Math.sqrt(u * u + v * v);

                if (norm<1e-4) {
                    // to prevent numerical instability in low field regions, keep moving
                    // in the same direction
                    u = u_last;
                    v = v_last;
                    norm = Math.sqrt(u * u + v * v);
                }
                
                stream_x += step * u / norm;
                stream_y += step * v / norm;

                x_steps.push(stream_x);
                y_steps.push(stream_y);

                // corrector step
                compute_field_electrostatic_direct_to_buffer(field, stream_x, stream_y, buffer);
                const u2 = buffer[0];
                const v2 = buffer[1];
                const norm2 = Math.sqrt(u2 * u2 + v2 * v2);

                if(norm2>1e-4) {
                    stream_x += step * (u2 / norm2 - u / norm)/2;
                    stream_y += step * (v2 / norm2 - v / norm)/2;
                }

                length_covered += Math.abs(step);
                
                u_last = u2;
                v_last = v2;
            }

            


            let landed_charge = getChargeFromPoint(charges, stream_x, stream_y);


            if (landed_charge !== null) {
                // the arrival at another charge needs to be registered so that we don't over-launch
                // field lines from that charge
                const angle = Math.atan2(-step*v_last, -step*u_last);
                debug_log("landed", landed_charge.id,"angle",angle);
                const arrived_at = landed_charge.departures;
                if(arrived_at.next_departure_index==arrived_at.num_departures) {
                    // Argh! We made a mistake. We've arrived at a charge that has no more field lines to launch.
                    // We need to backtrack and try again.
                    num_failed_launches++;
                    if(num_failed_launches>max_failed_launches) {
                        debug_log("Too many failed streamline launches. Will display the wrong number of fieldlines from at least one charge.");
                    } else {
                        debug_log("Failed streamline launch. Backtracking.");
                        departures.next_departure_index--;
                        departures.register_arrival(stream_angle, landed_charge, charge);
                        departures.num_departures++;
                        continue;
                    }

                } else {
                    arrived_at.register_arrival(angle, charge, landed_charge);
                }
            } else {
                debug_log("left the arena");
            }

            ctx.strokeStyle = 'black';
            ctx.lineWidth = 1.0;
            
            let i_start = -1;
            for(let i=0; i<x_steps.length; i++) {
                if(x_steps[i]>-step_size && y_steps[i]>-step_size && x_steps[i]<rect.width+step_size && y_steps[i]<rect.height+step_size) {
                    if(i_start == -1) {
                        i_start = i;
                        ctx.beginPath();
                        ctx.moveTo(x_steps[i], y_steps[i]);
                    } else {
                        ctx.lineTo(x_steps[i], y_steps[i]);
                    }
                } else if (i_start>=0) {
                    // we have just exited the screen, so draw an arrow in the middle of the segment just completed
                    ctx.stroke();
                    drawArrowOnFieldline(ctx, field, x_steps, y_steps, i_start, i);
                    i_start = -1 
                }
            }

            if(i_start>=0) {
                ctx.stroke();
                drawArrowOnFieldline(ctx, field, x_steps, y_steps, i_start, x_steps.length-1);
            }

                        
        }

    }       
    clearFieldlineDrawingInfo(charges); 
}

function drawArrowOnFieldline(ctx, field, x_steps, y_steps, i_start, i_end) {
    const index = Math.round((i_start + i_end)/2);
    let stream_x = x_steps[index];
    let stream_y = y_steps[index];

    const E = compute_field_electrostatic_direct(field, stream_x, stream_y);
    drawDirectionArrow(stream_x, stream_y, E.u, E.v, ctx);
}

function drawDirectionArrow(x_position, y_position, u, v, ctx) {
    ctx.save();
    ctx.translate(x_position, y_position);
    ctx.rotate(Math.atan2(v, u));
    ctx.beginPath();
    ctx.moveTo(-8, -5);
    ctx.lineTo(0, 0);
    ctx.lineTo(-8, 5);
    ctx.stroke();
    ctx.restore();
}


export function drawPotentialContours(charges, levels, ctx, color, show_direction = false) {
    let contours;
    if (show_direction) {
        let arrows;
        [contours, arrows] = generate_potential_contours_and_arrow_positions_at_levels(charges, levels);
        for (let arrow of arrows) {
            const [x, y] = arrow;
            const B = compute_field_magnetostatic_direct(charges, x, y);
            drawDirectionArrow(x, y, B.u, B.v, ctx);
        }
    } else {
        contours = generate_potential_contours_at_levels(charges, levels);
    }

    for (let contour of contours) {
        ctx.strokeStyle = color;
        ctx.lineWidth = 1.0;
        ctx.beginPath();
        for (let i = 0; i < contour.length; i++) {
            const [x, y] = contour[i];
            if (i === 0) {
                ctx.moveTo(x, y);
            } else {
                ctx.lineTo(x, y);
            }
        }
        ctx.stroke();
    }


}