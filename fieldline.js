import { compute_field_electrostatic_direct, generate_potential_contours_at_levels, generate_potential_contours_and_arrow_positions_at_levels, compute_field_magnetostatic_direct } from './maxwell/out/maxwell.js';
import { getChargeFromPoint, chargeSize } from './ui.js';

class MultiStreamDepartures {
    constructor(num_departures, starting_angle = 0) {
        this.num_departures = num_departures;
        this.next_departure_index = 0;
        this.starting_angle = starting_angle;
        this.arrival_range = {};
        this.arrivals = {}

    }

    get_new_departure_no_arrivals(departure_index_interleaved) {
        return this.starting_angle + (departure_index_interleaved+0.5) * 2*Math.PI / this.num_departures;
    }

    get_new_departure_with_arrivals(departure_index_interleaved) {
        let logical_angle = (departure_index_interleaved+0.5)/this.num_departures;


        // now the space available for actually launching streamlines is the space left after
        // excluding all arrival ranges. To complicate matters, the arrival_ranges may overlap, and
        // may also wrap around the 2 pi interval. So we need to find all the spaces left after 
        // excluding the arrival ranges, with these two complications in mind.
        let available_ranges = [{start: 0, end: 2*Math.PI}];
        let exclusion_ranges = []

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
                if (range.start<exclusion_start && range.end>exclusion_start) {
                    new_ranges.push({start: range.start, end: exclusion_start});
                }
                if (range.end>exclusion_end && range.start<exclusion_end) {
                    new_ranges.push({start: exclusion_end, end: range.end});
                }
                if ((exclusion_start<range.start && exclusion_end<range.start) 
                    || (exclusion_start>range.end && exclusion_end>range.end)) {
                    new_ranges.push(range);
                }
            }
            available_ranges = new_ranges;
        }
        console.log("prohibited:",exclusion_ranges);
        console.log("available:",available_ranges);
        console.log("offset:",this.starting_angle);

        
        
        // Now map the logical angle (0 to 1) into the actual angle space available:
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

        let index_within_range = departure_index_interleaved;
        let range_index = 0;
        while (index_within_range >= num_departure_per_range[range_index]) {
            index_within_range -= num_departure_per_range[range_index];
            range_index++;
        }
        

        console.log("num_departure_per_range:", num_departure_per_range);
        console.log("index_within_range:", index_within_range);
        console.log("range_index:", range_index);

        let range_min = available_ranges[range_index].start;
        let range_max = available_ranges[range_index].end;
        let angle = range_min + (range_max - range_min) * (index_within_range + 0.5) / num_departure_per_range[range_index];

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

    arrivals_to_range(arrivals) {
        let min_angle = Math.min(...arrivals);
        let max_angle = Math.max(...arrivals);
        let angle_range = max_angle - min_angle;


        // the range above may not be the best way to express it; for example if two lines
        // arrive, one at 6.1 and one at 0.1, we want to consider the arrival range to be
        // between 6.1 and 0.1 + 2 pi. Test for this condition and choose the best way
        // to express the interval
        for (let offset = 1; offset<4; offset++) {
            const arrivals_wrapped = arrivals.map(a => (a + offset * Math.PI/2)%(2*Math.PI) - offset * Math.PI/2);
            const min_angle_wrapped = Math.min(...arrivals_wrapped);
            const max_angle_wrapped = Math.max(...arrivals_wrapped);
            const offset_range = max_angle_wrapped - min_angle_wrapped;
            if (offset_range<angle_range) {
                min_angle = min_angle_wrapped;
                max_angle = max_angle_wrapped;
                angle_range = offset_range;
            }

        }

        console.log("min_angle:", min_angle, "max_angle:", max_angle, "angle_range:", angle_range);
        
        min_angle+=2*Math.PI;
        max_angle+=2*Math.PI;

        min_angle = min_angle%(2*Math.PI);
        max_angle = max_angle%(2*Math.PI);

        return [min_angle, max_angle];
    }

    register_arrival(angle, source) {
        angle = (angle + 2*Math.PI)%(2*Math.PI);
        console.log("registering arrival from", source, "at angle", angle);
        this.num_departures--;
        if (!(source in this.arrivals)) {
            this.arrivals[source] = [];
        }
        this.arrivals[source].push(angle);

        this.arrival_range[source] = this.arrivals_to_range(this.arrivals[source]);
        console.log("arrivals:", this.arrivals[source], "range:", this.arrival_range[source]);
        
    }
}

function fieldlineStartingAngles(charges) {
    // circle around the charge and measure the strength of the EM field
    const numCharges = charges.length;
    if (numCharges === 1) {
        charges[0].angle = 0;
    } else if (numCharges > 1) {
        for (let i = 0; i < numCharges; i++) {
            const currentCharge = charges[i];
            const nextCharge = charges[(i + 1) % numCharges];
            const dx = nextCharge.x - currentCharge.x;
            const dy = nextCharge.y - currentCharge.y;
            if (dx === 0 && dy === 0) {
                currentCharge.angle = 0;
            } else {
                currentCharge.angle = Math.atan2(dy, dx);
            }
        }
    }
}

export function drawElectrostaticFieldLines(charges, field, ctx, rect, chargeSize) {

    fieldlineStartingAngles(charges);
    let departures_all_charges = charges.map(charge => 
        ({charge: charge, departures: new MultiStreamDepartures(Math.abs(charge.charge) * 4,
                                                           charge.angle)}));

    // sort departures_all_charges by charge magnitude in ascending order
    departures_all_charges.sort((a, b) => Math.abs(a.charge.charge) - Math.abs(b.charge.charge));
    

    for (let {charge, departures} of departures_all_charges) {
        const x = charge.x;
        const y = charge.y;
        
        console.log("PROCESS CHARGE:",charge);

        let num_failed_launches = 0;
        const max_failed_launches = 20;

        while(true) {
            
            const stream_angle = departures.get_new_departure();
            if (stream_angle === null) break;

            console.log("Departure number", departures.next_departure_index, "of", departures.num_departures, "angle:", stream_angle);

            let stream_x = x + chargeSize * Math.cos(stream_angle);
            let stream_y = y + chargeSize * Math.sin(stream_angle);
            let length_covered = 0;

            
            

            

            const step = charge.charge>0?5.0:-5.0;

            const x_steps = [stream_x];
            const y_steps = [stream_y];

            let n_steps = 0;

            let u_last = 0;
            let v_last = 0;

            while((length_covered<20 || getChargeFromPoint(charges, stream_x, stream_y) === null)
                && (stream_x>0 && stream_y>0 && stream_x<rect.width && stream_y<rect.height)
                && n_steps<1000) {
                n_steps++;

                const E = compute_field_electrostatic_direct(field, stream_x, stream_y);
                let u = E.u;
                let v = E.v;
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
                const E2 = compute_field_electrostatic_direct(field, stream_x, stream_y);
                const u2 = E2.u;
                const v2 = E2.v;
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

            console.log("landed", landed_charge, "at", stream_x, stream_y)
            if (landed_charge !== null) {
                
                // the arrival at another charge needs to be registered so that we don't over-launch
                // field lines from that charge
                const angle = Math.atan2(stream_y - landed_charge.y, stream_x - landed_charge.x);
                const arrived_at = departures_all_charges.find(d => d.charge === landed_charge).departures;
                if(arrived_at.next_departure_index==arrived_at.num_departures) {
                    // Argh! We made a mistake. We've arrived at a charge that has no more field lines to launch.
                    // We need to backtrack and try again.
                    num_failed_launches++;
                    if(num_failed_launches>max_failed_launches) {
                        console.log("Too many failed streamline launches. Will display the wrong number of fieldlines from at least one charge.");
                    } else {
                        console.log("Failed streamline launch. Backtracking.");
                        departures.next_departure_index--;
                        departures.register_arrival(stream_angle, landed_charge.id);
                        departures.num_departures++;
                        continue;
                    }

                } else {
                    arrived_at.register_arrival(angle, charge.id);
                }
            }

            ctx.strokeStyle = 'black';
            ctx.lineWidth = 1.0;
            ctx.beginPath();
            ctx.moveTo(x_steps[0], y_steps[0]);
            for(let i=1; i<x_steps.length; i++)
                ctx.lineTo(x_steps[i], y_steps[i]);
            ctx.stroke();

            stream_x = x_steps[Math.round(x_steps.length/2)];
            stream_y = y_steps[Math.round(y_steps.length/2)];

            const E = compute_field_electrostatic_direct(field, stream_x, stream_y);
            drawDirectionArrow(stream_x, stream_y, E.u, E.v, ctx);            
        }

    }        
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