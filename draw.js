import { drawElectrostaticFieldLines, drawPotentialContours } from './fieldline.js';
import { compute_field_electrostatic_direct_to_buffer, compute_field_magnetostatic_direct_to_buffer, 
    compute_electric_field_dynamic_to_buffer, compute_field_electrostatic_per_charge_direct_to_buffer,
    compute_field_magnetostatic_per_charge_direct_to_buffer, 
    compute_one_force_electrostatic, compute_one_force_magnetostatic
} from './maxwell/out/maxwell.js';

export const chargeSize = 10;
const forceScaling = 0.1;

export function getChargeFromPoint(charges, x, y, allowRadius, addChargeSize=true, excludeCharge=null) {
    if(addChargeSize) {
        if (allowRadius == null) 
            allowRadius = chargeSize;
        else
            allowRadius += chargeSize;
    }

    for (let i = charges.length - 1; i >= 0; i--) {
        // go in reverse order so that the charge on top is selected first
        const charge = charges[i];
        const dx = x - charge.x;
        const dy = y - charge.y;
        if (Math.sqrt(dx * dx + dy * dy) < allowRadius && charge !== excludeCharge) {
            return charge;
        }
    }
    return null;
}

function drawChargesOrCurrents(ctx, charges, computeField, selectedCharge, forces) {
    const width = ctx.canvas.clientWidth;

    const charges_no_test_charges = charges.filter(charge => !charge.isTestCharge);

    if (forces !== null) {
        for (let i = 0; i < charges.length && i < forces.length; i++) {
            const charge = charges_no_test_charges[i];
            const force = forces[i];
            drawForce(ctx, charge, force);
        }
    }
    if (computeField === compute_field_magnetostatic_direct_to_buffer || computeField === compute_field_magnetostatic_per_charge_direct_to_buffer) {
        drawCurrents(ctx, charges, selectedCharge);
    } else {
        drawCharges(ctx, charges, selectedCharge);
    }
}

function drawForce(ctx, charge, force) {
    drawArrow(ctx, charge.x, charge.y, force.u*forceScaling*ctx.canvas.clientWidth, force.v*forceScaling*ctx.canvas.clientWidth, 'purple', 2, ctx.canvas.clientWidth/2, 20, false);
}

function drawTestChargeForces(ctx, charges, computeField, field) {
    charges.forEach(charge => {
        if(charge.isTestCharge) {
            let force;
            if (computeField === compute_field_electrostatic_direct_to_buffer) {
                force = compute_one_force_electrostatic(field, charge.x, charge.y, charge.charge);
            }
            else if (computeField === compute_field_magnetostatic_direct_to_buffer) {
                force = compute_one_force_magnetostatic(field, charge.x, charge.y, charge.charge);
            } else {
                force = {u: 0, v: 0};
            }

            drawForce(ctx, charge, force);
        }
    });
}

function drawCharges(ctx, charges, selectedCharge) {
    charges.forEach(charge => {
        if(charge.isTestCharge) {
            ctx.globalAlpha = 0.5;
        }
        
        ctx.beginPath();
        ctx.arc(charge.x, charge.y, chargeSize, 0, 2 * Math.PI, false);
        ctx.fillStyle = charge.charge > 0 ? 'red' : 'blue';
        ctx.fill();
        changeLineStyleIfSelected(ctx, charge, selectedCharge);
        ctx.stroke();
        if(charge.isTestCharge) {
            ctx.globalAlpha = 1;
        }
    });
}

function changeLineStyleIfSelected(ctx, charge, selectedCharge) {
    if (charge === selectedCharge) {
        ctx.strokeStyle = 'black';
        ctx.lineWidth = 2;
    } else {
        ctx.strokeStyle = 'grey';
        ctx.lineWidth = 1;
    }
}

function drawCurrents(ctx, charges, selectedCharge) {
    // Interpret the charges as currents. If current is positive, show as a 
    // circle with a dot in the centre. If current is negative, show as a circle with
    // a cross through it
    charges.forEach(charge => {
        ctx.save();
        if(charge.isTestCharge) { 
            ctx.strokeStlye = 'purple';
        } else {
            ctx.strokeStyle = 'black';
        }
        ctx.lineWidth = 1;

        ctx.beginPath();
        ctx.arc(charge.x, charge.y, chargeSize, 0, 2 * Math.PI, false);
        ctx.fillStyle = 'white';
        ctx.fill();
        // changeLineStyleIfSelected(ctx, charge, selectedCharge);
        ctx.stroke();

        if(charge.isTestCharge) { 
            ctx.fillStyle = 'purple';
        } else {
            ctx.fillStyle = 'black';
        }

        ctx.beginPath();
        if (charge.charge>0) {
            ctx.arc(charge.x, charge.y, chargeSize / 4, 0, 2 * Math.PI, false);
            
            ctx.fill();
            ctx.stroke();
        } else if (charge.charge < 0) {
            ctx.beginPath();
            ctx.moveTo(charge.x - chargeSize / 2, charge.y - chargeSize / 2);
            ctx.lineTo(charge.x + chargeSize / 2, charge.y + chargeSize / 2);
            ctx.stroke();
            ctx.beginPath();
            ctx.moveTo(charge.x + chargeSize / 2, charge.y - chargeSize / 2);
            ctx.lineTo(charge.x - chargeSize / 2, charge.y + chargeSize / 2);
            ctx.stroke();
        }
        ctx.restore();
    });
}

function generateVectors(computeField, rect, charges, field) {
    const vectors = [];
    const step = 20;
    var buffer;
    var n_per_point;
    var color = 'black';
    
    if(computeField === compute_field_electrostatic_per_charge_direct_to_buffer || computeField === compute_field_magnetostatic_per_charge_direct_to_buffer) {
        buffer = new Float64Array(2*charges.length);
        n_per_point = charges.length;
    } else {
        buffer = new Float64Array(2);
        n_per_point = 1;
    }

    for (let x = step; x < rect.width; x += step) {
        for (let y = step; y < rect.height; y += step) {
            // exclude the vector if it's within step distance from any charge
            if (charges.some(charge => {
                const dx = x - charge.x;
                const dy = y - charge.y;
                return Math.sqrt(dx * dx + dy * dy) < step;
            })) continue;
            computeField(field, x, y, buffer);
            for (let i = 0; i < n_per_point; i++) {
                if (n_per_point > 1) {
                    color = charges[i].charge > 0 ? 'red' : 'blue';
                }
                vectors.push({x, y, u: buffer[0+2*i], v: buffer[1+2*i], color});
            }
        }
    }
    
    return vectors;
}

function drawQuiverPlot(ctx, vectors) {
    // ctx.globalCompositeOperation = 'multiply';  // -- ideally this would be good for the 'show per charge' option but it seems slow
    vectors.forEach(({x, y, u, v, color}) => {
        drawArrow(ctx, x, y, u, v, color);
    });
    // ctx.globalCompositeOperation = 'source-over';
}


function drawArrow(ctx, x, y, u, v, color='black', linewidth=1, arrowLengthLimit=40, maxArrowHeadLength=8, centred=true) {
    let arrowLength = Math.sqrt(u * u + v * v);
    const angle = Math.atan2(v, u);

    let arrowHeadLength = arrowLength * 0.2;
    if(arrowHeadLength > maxArrowHeadLength) {
        arrowHeadLength = maxArrowHeadLength;
    }


    if(arrowLength > arrowLengthLimit) {
        u = u / arrowLength * arrowLengthLimit;
        v = v / arrowLength * arrowLengthLimit;
        arrowLength = arrowLengthLimit;
    }

    if(centred) {
        x-=u/2;
        y-=v/2;
    }

    ctx.lineWidth = linewidth;
    ctx.strokeStyle = color;

    ctx.beginPath();
    ctx.moveTo(x, y);
    ctx.lineTo(x + u, y + v);
    ctx.stroke();

    ctx.beginPath();
    ctx.fillStyle = color;
    ctx.moveTo(x + u, y + v);
    ctx.lineTo(x + u - arrowHeadLength * Math.cos(angle - Math.PI / 6), y + v - arrowHeadLength * Math.sin(angle - Math.PI / 6));
    ctx.lineTo(x + u - arrowHeadLength * Math.cos(angle + Math.PI / 6), y + v - arrowHeadLength * Math.sin(angle + Math.PI / 6));
    ctx.closePath();
    ctx.fill();
}

export function draw(ctx, rect, charges, field, fieldVisType, computeField, showPotential, selectedCharge, forces) {
    ctx.clearRect(0, 0, rect.width, rect.height);


    if (showPotential) {
        drawPotentialContours(field, [0], ctx, 'grey');
        drawPotentialContours(field, [250., 500., 750., 1000., 1250., 1500.,], ctx, 'blue');
        drawPotentialContours(field, [-250., -500.,-750., -1000., -1250., -1500.], ctx, 'red');
    }

    if (fieldVisType === 'quiver') {
        const vectors = generateVectors(computeField, rect, charges, field);
        drawQuiverPlot(ctx, vectors);
    } else if (fieldVisType === 'fieldline' && computeField === compute_field_electrostatic_direct_to_buffer) {
        const charges_no_test_charges = charges.filter(charge => !charge.isTestCharge);
        // The following algorithm only works when field lines start and end on charges, so perfect for the
        // electric case but not the magnetic case
        drawElectrostaticFieldLines(charges_no_test_charges, field, ctx, rect, chargeSize);
    } else if (fieldVisType === 'fieldline' && computeField === compute_field_magnetostatic_direct_to_buffer) {
        // Here we take cheeky advantage of the fact that the magnetostatic field lines are equivalent to
        // equipotential lines if we were solving an electrostatic problem. 
        const rangeValues = [];
        for (let i = 1.4; i <= 4.0; i += 0.4) {
            rangeValues.push(10**i);
            rangeValues.push(-(10**i));
        }
        rangeValues.push(0.0);
        drawPotentialContours(field, rangeValues, ctx, 'black', true);

    } else {
        console.error('Unknown field visualization type: ' + fieldVisType);
    }

    drawTestChargeForces(ctx, charges, computeField, field);
    drawChargesOrCurrents(ctx, charges, computeField, selectedCharge, forces);
}

    
