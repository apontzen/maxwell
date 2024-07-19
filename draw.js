import { drawElectrostaticFieldLines, drawPotentialContours } from './fieldline.js';
import { compute_field_electrostatic_direct, compute_field_magnetostatic_direct, 
    compute_electric_field_dynamic} from './maxwell/out/maxwell.js';


const maxArrowLength = 40;

export const chargeSize = 10;

export function getChargeFromPoint(charges, x, y, allowRadius) {
    if (allowRadius == null) 
        allowRadius = chargeSize;
    for (let i = charges.length - 1; i >= 0; i--) {
        // go in reverse order so that the charge on top is selected first
        const charge = charges[i];
        const dx = x - charge.x;
        const dy = y - charge.y;
        if (Math.sqrt(dx * dx + dy * dy) < allowRadius) {
            return charge;
        }
    }
    return null;
}

function drawChargesOrCurrents(ctx, charges, computeField, selectedCharge) {
    if (computeField === compute_field_magnetostatic_direct) {
        drawCurrents(ctx, charges, selectedCharge);
    } else {
        drawCharges(ctx, charges, selectedCharge);
    }
}

function drawCharges(ctx, charges, selectedCharge) {
    charges.forEach(charge => {
        
        ctx.beginPath();
        ctx.arc(charge.x, charge.y, chargeSize, 0, 2 * Math.PI, false);
        ctx.fillStyle = charge.charge > 0 ? 'red' : 'blue';
        ctx.fill();
        changeLineStyleIfSelected(charge, selectedCharge);
        ctx.stroke();
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
        ctx.beginPath();
        ctx.arc(charge.x, charge.y, chargeSize, 0, 2 * Math.PI, false);
        ctx.fillStyle = 'white';
        ctx.fill();
        changeLineStyleIfSelected(charge, selectedCharge);
        ctx.stroke();
        ctx.beginPath();
        if (charge.charge>0) {
            ctx.arc(charge.x, charge.y, chargeSize / 4, 0, 2 * Math.PI, false);
            ctx.fillStyle = 'black';
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
    });
}

function generateVectors(computeField, rect, charges, field) {
    const vectors = [];
    const step = 20;

    for (let x = step; x < rect.width; x += step) {
        for (let y = step; y < rect.height; y += step) {
            // exclude the vector if it's within step distance from any charge
            if (charges.some(charge => {
                const dx = x - charge.x;
                const dy = y - charge.y;
                return Math.sqrt(dx * dx + dy * dy) < step;
            })) continue;
            const vector = computeField(field, x, y);
            vectors.push({x, y, u: vector.u, v: vector.v});
        }
    }
    
    return vectors;
}

function drawQuiverPlot(ctx, vectors) {
    vectors.forEach(({x, y, u, v}) => {
        drawArrow(ctx, x, y, u, v);
    });
}


function drawArrow(ctx, x, y, u, v) {
    let arrowLength = Math.sqrt(u * u + v * v);
    const angle = Math.atan2(v, u);

    if(arrowLength > maxArrowLength) {
        u = u / arrowLength * maxArrowLength;
        v = v / arrowLength * maxArrowLength;
        arrowLength = maxArrowLength;
    }

    x-=u/2;
    y-=v/2;

    ctx.lineWidth = 1;
    ctx.strokeStyle = 'black';

    ctx.beginPath();
    ctx.moveTo(x, y);
    ctx.lineTo(x + u, y + v);
    ctx.stroke();

    ctx.beginPath();
    ctx.fillStyle = 'black';
    ctx.moveTo(x + u, y + v);
    ctx.lineTo(x + u - arrowLength * 0.2 * Math.cos(angle - Math.PI / 6), y + v - arrowLength * 0.2 * Math.sin(angle - Math.PI / 6));
    ctx.lineTo(x + u - arrowLength * 0.2 * Math.cos(angle + Math.PI / 6), y + v - arrowLength * 0.2 * Math.sin(angle + Math.PI / 6));
    ctx.closePath();
    ctx.fill();
}

export function draw(ctx, rect, charges, field, fieldVisType, computeField, showPotential, selectedCharge) {
    ctx.clearRect(0, 0, rect.width, rect.height);


    if (showPotential) {
        drawPotentialContours(field, [0], ctx, 'grey');
        drawPotentialContours(field, [125., 250., 500., 1000.], ctx, 'blue');
        drawPotentialContours(field, [-125., -250., -500., -1000.], ctx, 'red');
    }

    if (fieldVisType === 'quiver') {
        const vectors = generateVectors(computeField, rect, charges, field);
        drawQuiverPlot(ctx, vectors);
    } else if (fieldVisType === 'fieldline' && computeField === compute_field_electrostatic_direct) {
        // The following algorithm only works when field lines start and end on charges, so perfect for the
        // electric case but not the magnetic case
        drawElectrostaticFieldLines(charges, field, ctx, rect, chargeSize);
    } else if (fieldVisType === 'fieldline' && computeField === compute_field_magnetostatic_direct) {
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

    drawChargesOrCurrents(ctx, charges, computeField, selectedCharge);
}

    
