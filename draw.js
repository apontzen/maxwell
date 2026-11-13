import { drawPotentialContours, drawPotentialContours2D, drawElectricFieldLines2D,
    drawMagneticFieldLines2D } from './fieldline.js';

export const chargeSize = 10;
const forceScaling = 1.0/200**3;

export function getChargeFromPoint(charges, x, y, allowRadius, addChargeSize=true, excludeCharge=null, dipoleMode=false) {
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

    if (dipoleMode) {
        for (let i = charges.length - 1; i >= 0; i--) {
            const charge = charges[i];
            if (charge.dipoleWith !== undefined && charge.dipoleWith > charge.id) {
                const otherCharge = charges.find(c => c.id === charge.dipoleWith);
                // find the closest point to (x,y) on the line between the two charges
                const dx = otherCharge.x - charge.x;
                const dy = otherCharge.y - charge.y;
                const t = ((x - charge.x) * dx + (y - charge.y) * dy) / (dx * dx + dy * dy);
                const closestX = charge.x + t * dx;
                const closestY = charge.y + t * dy;
                const distSquared = (x - closestX) ** 2 + (y - closestY) ** 2;
                if (t >= 0 && t <= 1 && distSquared < allowRadius**2) {
                    return [charge, otherCharge];
                }

            }
        }
    }
    return null;
}

function drawChargesOrCurrents(ctx, charges, magnetic, selectedCharge, forces, dipoleMode) {
    const width = ctx.canvas.clientWidth;

    const charges_no_test_charges = charges.filter(charge => !charge.isTestCharge);

    if (forces !== null) {
        for (let i = 0; i < charges.length && i < forces.length; i++) {
            const charge = charges_no_test_charges[i];
            const force = forces[i];
            drawForce(ctx, charge.x, charge.y, force);
        }
    }
    if(dipoleMode) {
        charges.forEach(charge => {
            if(charge.dipoleWith !== undefined) {
                if(charge.dipoleWith > charge.id) {
                    const otherCharge = charges.find(c => c.id === charge.dipoleWith);
                    ctx.beginPath();
                    ctx.moveTo(charge.x, charge.y);
                    ctx.lineTo(otherCharge.x, otherCharge.y);
                    if(charge.isTestCharge) {
                        ctx.strokeStyle = 'purple';
                    } else {
                        ctx.strokeStyle = 'grey';
                    }
                    ctx.lineWidth = 2;
                    ctx.stroke();
                }
            }
        });
    }

    if (magnetic) {
        drawCurrents(ctx, charges, selectedCharge);
    } else {
        drawCharges(ctx, charges, selectedCharge);
    }
}

function drawForce(ctx, x, y, force) {
    // We scale the forces by the canvas size cubed. This might look weird at first, but it's because natively the forces are calculated as inverse square laws
    // in pixel space, so to make them scale down we need a canvas size squared factor; then another factor of the canvas size to scale them to the same
    // length relative to the canvas.
    drawArrow(ctx, x, y, force.u*forceScaling*ctx.canvas.clientWidth**3, force.v*forceScaling*ctx.canvas.clientWidth**3, 'purple', 2, ctx.canvas.clientWidth/2, 20, false, 5);
}

function drawTorque(ctx, x_origin, y_origin, x_start, y_start, num_degrees) {
    // Draw a circle segment to represent the torque, starting at x_start, y_start and going num_degrees clockwise around x_origin, y_origin
    const radius = Math.sqrt((x_start - x_origin) ** 2 + (y_start - y_origin) ** 2);
    const angle_start = Math.atan2(y_start - y_origin, x_start - x_origin);
    const angle_end = angle_start + num_degrees * Math.PI / 180;

    ctx.beginPath();
    if(angle_start<angle_end) {
        ctx.arc(x_origin, y_origin, radius, angle_start, angle_end);
    } else {
        ctx.arc(x_origin, y_origin, radius, angle_end, angle_start);
    }
    ctx.strokeStyle = 'purple';
    ctx.lineWidth = 2;
    ctx.stroke();



    let arrow_angle = angle_end + Math.PI/2;
    let arrow_offset = 2.0/radius;
    if (num_degrees < 0) {
        arrow_offset = -arrow_offset;
        arrow_angle = angle_end - Math.PI/2;
    }

    drawArrowHead(ctx, x_origin + radius * Math.cos(angle_end+arrow_offset), y_origin + radius * Math.sin(angle_end+arrow_offset), 10, arrow_angle, 'purple');
    
    



}

function logScale(x, minval, maxval) {
    let abs_x = Math.abs(x);
    let sign_x = Math.sign(x);
    if (abs_x < minval) {
        return x;
    }
    let log_x = Math.log(abs_x);
    let log_minval = Math.log(minval);
    let log_maxval = Math.log(maxval);
    let log_range = log_maxval - log_minval;
    let log_x_scaled = (log_x - log_minval) / log_range;
    return sign_x * maxval * log_x_scaled;
}


function drawTestChargeForces(ctx, charges, renderModel, dipoleMode) {
    charges.forEach(charge => {
        if(charge.isTestCharge) {
            let force;
            // In per-charge mode no forces are shown (as before the field-model refactor)
            const fieldComputationFunction = renderModel.perCharge ? () => ({u: 0, v: 0})
                : (x, y, q) => renderModel.forceOn(x, y, q);

            force = fieldComputationFunction(charge.x, charge.y, charge.charge);
            let x = charge.x;
            let y = charge.y;

            if(dipoleMode && charge.dipoleWith !== undefined) {
                if (charge.dipoleWith > charge.id) {
                    const otherCharge = charges.find(c => c.id === charge.dipoleWith);
                    x = (x + otherCharge.x) / 2;
                    y = (y + otherCharge.y) / 2;
                    const otherForce = fieldComputationFunction(otherCharge.x, otherCharge.y, otherCharge.charge);
                    let torque_around_com = force.u * (otherCharge.y - y) - force.v * (otherCharge.x - x);
                    torque_around_com += otherForce.u * (charge.y - y) - otherForce.v * (charge.x - x);

                    const distanceBetweenCharges = Math.sqrt((otherCharge.x - charge.x)**2 + (otherCharge.y - charge.y)**2);
                    torque_around_com *= 5.0e-2*ctx.canvas.clientWidth/(distanceBetweenCharges); // so that length of torque line, not angle, represents torque
                
                    
                    if (torque_around_com > 160) {
                        torque_around_com = 160;
                    } else if (torque_around_com < -160) {
                        torque_around_com = -160;
                    }
        

                    drawTorque(ctx, x, y, otherCharge.x, otherCharge.y, torque_around_com);
                    drawTorque(ctx, x, y, charge.x, charge.y, torque_around_com);
        
                    force.u += otherForce.u;
                    force.v += otherForce.v;
                } else {
                    return;
                }
            }

            drawForce(ctx, x, y, force);
        }
    });
}

function drawCharges(ctx, charges, selectedCharge) {
    charges.forEach(charge => {        
        ctx.beginPath();
        ctx.arc(charge.x, charge.y, chargeSize, 0, 2 * Math.PI, false);
        if(charge.isTestCharge) {
            ctx.fillStyle = charge.charge >0 ? 'lightcoral' : 'lightblue';
        } else {
            ctx.fillStyle = charge.charge > 0 ? 'red' : 'blue';
        }
        ctx.fill();
        changeLineStyleIfSelected(ctx, charge, selectedCharge);
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
        ctx.save();
        if(charge.isTestCharge) { 
            ctx.strokeStyle = 'purple';
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

function generateVectors(renderModel, rect, charges) {
    const vectors = [];
    const step = 20;
    var buffer;
    var n_per_point;
    var color = 'black';
    
    const field = renderModel.field;
    let computeField = renderModel.fieldToBuffer;
    if(renderModel.perCharge) {
        computeField = renderModel.fieldPerChargeToBuffer;
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

function drawArrowHead(ctx, x, y, arrowHeadLength, angle, color) {
    ctx.beginPath();
    ctx.fillStyle = color;
    ctx.moveTo(x, y);
    ctx.lineTo(x - arrowHeadLength * Math.cos(angle - Math.PI / 6), y - arrowHeadLength * Math.sin(angle - Math.PI / 6));
    ctx.lineTo(x - arrowHeadLength * Math.cos(angle + Math.PI / 6), y - arrowHeadLength * Math.sin(angle + Math.PI / 6));
    ctx.closePath();
    ctx.fill();
}


export function drawArrow(ctx, x, y, u, v, color='black', linewidth=1, arrowLengthLimit=40, maxArrowHeadLength=8, centred=true, minArrowHeadLength=0) {
    let arrowLength = Math.sqrt(u * u + v * v);
    const angle = Math.atan2(v, u);

    let arrowHeadLength = arrowLength * 0.2;
    if(arrowHeadLength > maxArrowHeadLength) {
        arrowHeadLength = maxArrowHeadLength;
    }
    if(arrowHeadLength < minArrowHeadLength) {
        arrowHeadLength = minArrowHeadLength;
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

    let arrowHeadX = x + u + 0.5*arrowHeadLength*u/arrowLength;
    let arrowHeadY = y + v + 0.5*arrowHeadLength*v/arrowLength;

    if (arrowHeadX<0) {
        let t = arrowHeadX/u;
        arrowHeadX = 0;
        arrowHeadY = arrowHeadY - t*v;
    }
    if (arrowHeadY<0) {
        let t = arrowHeadY/v;
        arrowHeadY = 0;
        arrowHeadX = arrowHeadX - t*u;
    }
    if (arrowHeadX>ctx.canvas.clientWidth) {
        let t = (arrowHeadX-ctx.canvas.clientWidth)/u;
        arrowHeadX = ctx.canvas.clientWidth;
        arrowHeadY = arrowHeadY - t*v;
    }
    if (arrowHeadY>ctx.canvas.clientHeight) {
        let t = (arrowHeadY-ctx.canvas.clientHeight)/v;
        arrowHeadY = ctx.canvas.clientHeight;
        arrowHeadX = arrowHeadX - t*u;
    }

    drawArrowHead(ctx, arrowHeadX, arrowHeadY, arrowHeadLength, angle, color);
}

/**
 * Draw the whole scene.
 *
 * @param renderModel  from solvers.js makeRenderModel(): the (effective) field model
 * @param options      {plotType: 'quiver'|'fieldline', showPotential, selectedCharge, forces,
 *                      dipoleMode, extraDrawInfo: {gaussianSurfacePoints}}
 */
export function draw(ctx, rect, charges, renderModel, { plotType, showPotential, selectedCharge, forces, dipoleMode, extraDrawInfo }) {
    ctx.clearRect(0, 0, rect.width, rect.height);

    const gaussianSurfacePoints = extraDrawInfo !== undefined ? extraDrawInfo['gaussianSurfacePoints'] : [];
    
    if (gaussianSurfacePoints.length > 0) {
        ctx.save();
        ctx.beginPath();
        ctx.strokeStyle = 'purple';
        ctx.fillStyle = 'lightgrey';
        ctx.lineWidth = 2;
        ctx.setLineDash([5, 5]);
        ctx.moveTo(gaussianSurfacePoints[0].x, gaussianSurfacePoints[0].y);
        for (let i = 1; i < gaussianSurfacePoints.length; i++) {
            ctx.lineTo(gaussianSurfacePoints[i].x, gaussianSurfacePoints[i].y);
        }
        ctx.closePath();
        

        ctx.stroke();
        ctx.fill();
        
        
        ctx.restore();
    }


    const model2D = renderModel.model2D;

    if (showPotential && !renderModel.isMagnetic && !renderModel.perCharge && renderModel.info.potential) {
        if (model2D) {
            drawPotentialContours2D(ctx, model2D, rect);
        } else {
            const field = renderModel.field;
            drawPotentialContours(field, [0], ctx, 'grey');
            drawPotentialContours(field, [250., 500., 750., 1000., 1250., 1500.,], ctx, 'blue');
            drawPotentialContours(field, [-250., -500.,-750., -1000., -1250., -1500.], ctx, 'red');
        }
    }

    if (plotType === 'quiver') {
        const vectors = generateVectors(renderModel, rect, charges);
        drawQuiverPlot(ctx, vectors);
    } else if (plotType === 'fieldline' && model2D && !model2D.isMagnetic) {
        drawElectricFieldLines2D(ctx, model2D, rect, gaussianSurfacePoints);
    } else if (plotType === 'fieldline' && model2D && model2D.isMagnetic) {
        drawMagneticFieldLines2D(ctx, model2D, rect, gaussianSurfacePoints);
    } else {
        console.error(`Cannot draw ${plotType} for solver ${renderModel.solver}`);
    }

    drawTestChargeForces(ctx, charges, renderModel, dipoleMode);
    drawChargesOrCurrents(ctx, charges, renderModel.isMagnetic, selectedCharge, forces, dipoleMode);
}

    
