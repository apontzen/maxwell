import { compute_field_electrostatic_direct } from './maxwell/out/maxwell.js';
import { getChargeFromPoint, chargeSize } from './ui.js';

class StreamDepartures {
    constructor(num_departures, starting_angle = 0) {
        this.num_departures = num_departures;
        this.next_departure_index = 0;
        this.starting_angle = starting_angle;
        this.ending_angle = starting_angle + 2*Math.PI;
        this.restricted_angle = false;
        this.arrivals = new Array();
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

        if (this.restricted_angle) {
            angle = this.starting_angle + (departure_index_interleaved+1) * (this.ending_angle - this.starting_angle) / (this.num_departures + 1);
        } else {
            angle = this.starting_angle + (departure_index_interleaved+0.5) * (this.ending_angle - this.starting_angle) / (this.num_departures );
        }
        this.next_departure_index++;
        return angle%(2*Math.PI);
    }

    

    register_arrival(angle) {
        this.arrivals.push((angle+2*Math.PI)%(2*Math.PI));
        if (this.num_departures > 0)
            this.num_departures--;

        this.restricted_angle = this.arrivals.length > 1;

        let min_angle = Math.min(...this.arrivals);
        let max_angle = Math.max(...this.arrivals);
        let angle_range = max_angle - min_angle;


        // the range above may not be the best way to express it; for example if two lines
        // arrive, one at 6.1 and one at 0.1, we want to consider the arrival range to be
        // between 6.1 and 0.1 + 2 pi. Test for this condition and choose the best way
        // to express the interval
        for (let offset = 1; offset<4; offset++) {
            const arrivals_wrapped = this.arrivals.map(a => (a + offset * Math.PI/2)%(2*Math.PI) - offset * Math.PI/2);
            const min_angle_wrapped = Math.min(...arrivals_wrapped);
            const max_angle_wrapped = Math.max(...arrivals_wrapped);
            const offset_range = max_angle_wrapped - min_angle_wrapped;
            if (offset_range<angle_range) {
                min_angle = min_angle_wrapped;
                max_angle = max_angle_wrapped;
                angle_range = offset_range;
            }

        }

        this.starting_angle = max_angle;
        this.ending_angle = min_angle+2*Math.PI;

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

export function drawfieldlinePlot(charges, field, ctx, rect, chargeSize) {

    fieldlineStartingAngles(charges);
    let departures_all_charges = charges.map(charge => 
        ({charge: charge, departures: new StreamDepartures(Math.abs(charge.charge) * 4,
                                                           charge.angle)}));

    // sort departures_all_charges by charge magnitude in ascending order
    departures_all_charges.sort((a, b) => Math.abs(a.charge.charge) - Math.abs(b.charge.charge));
    

    for (let {charge, departures} of departures_all_charges) {
        const x = charge.x;
        const y = charge.y;

        let num_failed_launches = 0;
        const max_failed_launches = 20;

        while(true) {
            const stream_angle = departures.get_new_departure();
            if (stream_angle === null) break;

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
                let u = E[0];
                let v = E[1];
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
                const u2 = E2[0];
                const v2 = E2[1];
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
                        departures.register_arrival(stream_angle);
                        departures.num_departures++;
                        continue;
                    }

                } else {
                    arrived_at.register_arrival(angle);
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
            const u = E[0];
            const v = E[1];
            ctx.save();
            ctx.translate(stream_x, stream_y);
            ctx.rotate(Math.atan2(v, u));
            ctx.beginPath();
            ctx.moveTo(-8, -5);
            ctx.lineTo(0, 0);
            ctx.lineTo(-8, 5);
            //ctx.closePath();
            //ctx.fillStyle = 'black';
            //ctx.fill();
            ctx.stroke();
            ctx.restore();
            
        }

    }        
}