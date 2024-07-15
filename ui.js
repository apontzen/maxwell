import init, { compute_field_electrostatic_direct, compute_field_magnetostatic_direct, 
    compute_electric_field_dynamic, init_panic_hook, FieldConfiguration } from './maxwell/out/maxwell.js';

import { drawElectrostaticFieldLines, drawPotentialContours } from './fieldline.js';

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

export async function main() {
    await init();

    init_panic_hook();


    const canvas = document.getElementById('vectorFieldCanvas');
    const ctx = canvas.getContext('2d');
    const addPositiveChargeBtn = document.getElementById('addPositiveCharge');
    const addNegativeChargeBtn = document.getElementById('addNegativeCharge');
    // const fieldStrength = document.getElementById('fieldStrength');

    
    const maxArrowLength = 40;

    let charges = [];
    let draggingCharge = null;
    let draggingOffsetX = 0;
    let draggingOffsetY = 0;

    let nextChargeStrength = 1;

    let cic_resolution = 128;

    const dpr = window.devicePixelRatio || 1;

    // make the DOM element smaller if it's exceeding the width of the device
    canvas.style.width = '90%';
    canvas.style.height = '90%';
    canvas.style.maxWidth = '800px';
    canvas.style.maxHeight = '600px';
    
    // Store the original CSS dimensions.
    const rect = canvas.getBoundingClientRect();

    // Give the canvas pixel dimensions of their CSS
    // size * the device pixel ratio.
    canvas.width = rect.width * dpr;
    canvas.height = rect.height * dpr;

    // Scale all drawing operations by the dpr, so you
    // don't have to worry about the difference.
    ctx.scale(dpr, dpr);

    // Set the CSS dimensions to the original values.
    canvas.style.width = `${rect.width}px`;
    canvas.style.height = `${rect.height}px`;

    let computeField = compute_field_electrostatic_direct;
    let field = null;

    let dynamic = false;
    
    let plotType = 'quiver';

    function updateSolverType() {
        let allowPotential = false;

        // get value of solver
        const solverType = document.getElementById('solver').value;
        
        dynamic = solverType === 'dynamic';

        plotType = 'quiver';

        if (solverType.endsWith('_fieldline')) {
            plotType = 'fieldline';
        }

        if (solverType === 'electrostatic_direct' || solverType === 'electrostatic_direct_fieldline') {
            updateChargeOrCurrentLabel('Charge');
            computeField = compute_field_electrostatic_direct;
            allowPotential = true;
        } else if (solverType === 'electrostatic_fourier' || solverType === 'dynamic') {
            updateChargeOrCurrentLabel('Charge');
            computeField = compute_electric_field_dynamic;
        } else if (solverType === 'magnetostatic_direct' || solverType === 'magnetostatic_direct_fieldline') { 
            updateChargeOrCurrentLabel('Current');
            computeField = compute_field_magnetostatic_direct;
        } else {
            console.error('Unknown solver type');
        }
        field = new FieldConfiguration(rect.width, rect.height, cic_resolution, cic_resolution);

        if(allowPotential) {
            document.getElementById('show-potential-control').style.display = 'inline';
        } else {
            document.getElementById('show-potential-control').style.display = 'none';
        }
    }

    updateSolverType();

    let animation_request_id = null;
    let last_time = null;

    document.getElementById('solver').addEventListener('change', () => {
        if(animation_request_id!==null)
            window.cancelAnimationFrame(animation_request_id);
        animation_request_id = null;
        last_time = null;
        updateSolverType();
        drawVectorField();
        if(dynamic)
            animation_request_id = window.requestAnimationFrame(tickField);
    });

    document.getElementById('potential').addEventListener('change', () => {
        drawVectorField();
    });

    function stateToJson() {
        let charges_normalized_coordinates = charges.map(charge => {
            const normalizedX = charge.x / rect.width;
            const normalizedY = charge.y / rect.height;
            return { ...charge, x: normalizedX, y: normalizedY };
        });

        return JSON.stringify({
            charges: charges_normalized_coordinates,
            solver: document.getElementById('solver').value,
            show_potential: document.getElementById('potential').checked
        });
    }

    function jsonToState(str) {
        const state = JSON.parse(str);
        if (state) {
            charges = state.charges.map(charge => {
                const scaledX = charge.x * rect.width;
                const scaledY = charge.y * rect.height;
                return {...charge, x: scaledX, y: scaledY};
            });

            document.getElementById('solver').value = state.solver;
            document.getElementById('potential').checked = state.show_potential;
            updateSolverType();
            return true;
        } else {
            return false;
        }
    }

    function updateDisplayedJson() {
        document.getElementById('jsonDescription').value = stateToJson();
    }

    document.getElementById('loadJson').addEventListener('click', () => {
        const json = document.getElementById('jsonDescription').value;
        jsonToState(json);
        drawVectorField();
    });

    document.getElementById('copyJson').addEventListener('click', () => {
        const json = document.getElementById('jsonDescription').value;
        navigator.clipboard.writeText(json);
    });

    document.getElementById('pasteJson').addEventListener('click', async () => {
        const json = await navigator.clipboard.readText();
        console.log("paste", json);
        if(jsonToState(json)) {
            document.getElementById('jsonDescription').value = json;
            drawVectorField();
        }
    });

    // Load the state from local storage if it exists
    const savedState = localStorage.getItem('maxwell_state');
    if (savedState) {
        jsonToState(savedState);
    }

    // Save the state to local storage whenever it changes
    function saveState() {
        const state = stateToJson();
        localStorage.setItem('maxwell_state', state);
    }


    function tickField(time_now) {
        animation_request_id = null;

        let dt = 0;
        if (last_time === null) {
            last_time = time_now;
        } else {
            dt = time_now - last_time;
        }

        dt /= 5; // convert from ms to internal time units

        // Try to get forward in time by dt, but abandon if it takes too long (more than 20ms)
        const performance_time_start = performance.now();
        const max_time_ms = 20.0;

        while (dt > 0.0 && performance.now() - performance_time_start < max_time_ms) {
            const step = Math.min(0.5, dt);
            field.tick(step);
            dt -= step;
        }
 
        last_time = time_now;
        drawVectorField(false);
 

        // only request the next frame if we're still solving dynamically
        if(dynamic)
            animation_request_id = window.requestAnimationFrame(tickField);

    }

    function drawVectorField(user_update = true) {
        if(user_update) {
            updateDisplayedJson();
            saveState();
        }

        ctx.clearRect(0, 0, canvas.width, canvas.height);

        field.set_charges(charges);
        if(!dynamic) {
            field.reset_fields();
            // if dynamic, the fields are being updated in real time and we don't want to reset them
        }

        const x0 = 200.0;
        const y0 = 200.0;
        if (charges.length > 0 && computeField === compute_field_electrostatic_direct && document.getElementById('potential').checked) {
            drawPotentialContours(field, [0], ctx, 'grey');

            drawPotentialContours(field, [125., 250., 500., 1000.], ctx, 'blue');
            
            drawPotentialContours(field, [-125., -250., -500., -1000.], ctx, 'red');
        }

        if (plotType === 'quiver') {
            const vectors = generateVectors();
            drawQuiverPlot(vectors);
        } else if (plotType === 'fieldline') {
            if(computeField === compute_field_electrostatic_direct) {
                // The following algorithm only works when field lines start and end on charges, so perfect for the
                // electric case but not the magnetic case
                drawElectrostaticFieldLines(charges, field, ctx, rect, chargeSize);
            } else if (computeField === compute_field_magnetostatic_direct) {
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
                console.error('Fieldlines not supported for this solver');
            }
        }

        drawChargesOrCurrents();

    }

    function updateChargeOrCurrentLabel(charge_or_current) {
        // find all spans with class charge_or_current and update their contents to the provided string
        const chargeOrCurrentSpans = document.querySelectorAll(`.charge_or_current`);
        chargeOrCurrentSpans.forEach(span => {
            span.textContent = charge_or_current;
        });
    }

    function generateVectors() {
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

    function drawQuiverPlot(vectors) {
        vectors.forEach(({x, y, u, v}) => {
            drawArrow(x, y, u, v);
        });
    }


    function drawArrow(x, y, u, v) {
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

    const chargeInput = document.getElementById('charge');
    const chargeValue = document.getElementById('chargeValue');
    let selectedCharge = null;

    chargeInput.addEventListener('input', () => {
        selectedCharge.charge = parseInt(chargeInput.value);
        chargeValue.textContent = chargeInput.value;
        drawVectorField();
    });

    function deleteCharge(charge) {
        charges = charges.filter(c => c !== charge);
        field.set_charges(charges);
    }

    const deleteChargeBtn = document.getElementById('deleteCharge');
    deleteChargeBtn.addEventListener('click', () => {
        deleteCharge(selectedCharge);
        deselectCharge();
        drawVectorField();
    });


    function selectCharge(charge) {

        chargeInput.value = charge.charge;
        chargeValue.textContent = chargeInput.value;
        selectedCharge = charge;

        drawVectorField();

        const chargeProperties = document.querySelector('.charge-properties');
        chargeProperties.style.display = 'block';

        chargeProperties.style.display = 'block';
        chargeProperties.style.position = 'absolute';
        const canvasRect = canvas.getBoundingClientRect();
        const pageOffsetTop = window.scrollY || document.documentElement.scrollTop || document.body.scrollTop || 0;
        const pageOffsetLeft = window.scrollX || document.documentElement.scrollLeft || document.body.scrollLeft || 0;
        const canvasRectTop = canvasRect.top + pageOffsetTop;
        const canvasRectLeft = canvasRect.left + pageOffsetLeft;
        chargeProperties.style.top = `${canvasRectTop + charge.y - chargeProperties.offsetHeight / 2}px`;

        // Check if charge is on the right side of the screen
        if (charge.x > rect.width / 2) {
            chargeProperties.style.left = `${canvasRectLeft + charge.x - chargeProperties.offsetWidth - 20}px`;
            chargeProperties.classList.add('point-right');
            chargeProperties.classList.remove('point-left');
        } else {
            chargeProperties.style.left = `${canvasRectLeft + charge.x + 20}px`;
            chargeProperties.classList.add('point-left');
            chargeProperties.classList.remove('point-right');
        }
        // if the chargeProperties div is off the screen, move it back on
        if (chargeProperties.offsetLeft < 0) {
            chargeProperties.style.left = '0';
        }
        if (chargeProperties.offsetLeft + chargeProperties.offsetWidth > window.innerWidth) {
            chargeProperties.style.left = `${window.innerWidth - chargeProperties.offsetWidth}px`;
        }
        chargeProperties.style.zIndex = '1';

        
    }

    function deselectCharge() {
        const chargeProperties = document.querySelector('.charge-properties');
        chargeProperties.style.display = 'none';
        if (selectedCharge !== null) {
            nextChargeStrength = selectedCharge.charge;         
        }

        selectedCharge = null;
        drawVectorField();
    }

    function drawChargesOrCurrents() {
        if (computeField === compute_field_magnetostatic_direct) {
            drawCurrents();
        } else {
            drawCharges();
        }
    }

    function drawCharges() {
        charges.forEach(charge => {
            
            ctx.beginPath();
            ctx.arc(charge.x, charge.y, chargeSize, 0, 2 * Math.PI, false);
            ctx.fillStyle = charge.charge > 0 ? 'red' : 'blue';
            ctx.fill();
            changeLineStyleIfSelected(charge);
            ctx.stroke();
        });
    }

    function changeLineStyleIfSelected(charge) {
        if (charge === selectedCharge) {
            ctx.strokeStyle = 'black';
            ctx.lineWidth = 2;
        } else {
            ctx.strokeStyle = 'grey';
            ctx.lineWidth = 1;
        }
    }

    function drawCurrents() {
        // Interpret the charges as currents. If current is positive, show as a 
        // circle with a dot in the centre. If current is negative, show as a circle with
        // a cross through it
        charges.forEach(charge => {
            ctx.beginPath();
            ctx.arc(charge.x, charge.y, chargeSize, 0, 2 * Math.PI, false);
            ctx.fillStyle = 'white';
            ctx.fill();
            changeLineStyleIfSelected(charge);
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

    let charge_id = 0;

    function addCharge(x, y, charge) {
        charges.push({x, y, charge, id: charge_id});
        charge_id++;
        deselectCharge();
        field.set_charges(charges);
        drawVectorField();
    }

    function outputNumpyArray() {
        // sample the field along the x axis, with y set equal to charges[0].y, and
        // output to a string that can be parsed by numpy for debug purposes
        const y = charges[0].y;
        const x_values = [];
        const u_values = [];
        const v_values = [];
        for (let x = 0; x < rect.width; x += 10) {
            const E = computeField(field, x, y);
            x_values.push(x);
            u_values.push(E[0]);
            v_values.push(E[1]);
        }

        const x_ar = ('x = np.array(' + JSON.stringify(x_values) + ');');
        const u_ar = ('u = np.array(' + JSON.stringify(u_values) + ');');
        const v_ar = ('v = np.array(' + JSON.stringify(v_values) + ')');

        // set the contents of dom element id numpy-output to x_ar+u_ar+v_ar
        
        console.log(x_ar);
        console.log(u_ar);
        console.log(v_ar);

        navigator.clipboard.writeText(x_ar + u_ar + v_ar);
    }


    let originalChargeX = 0;
    let originalChargeY = 0;


    function coordinatesFromMouseOrTouch(event) {
        if (event.touches) {
            const touch = event.touches[0];
            return { offsetX: touch.pageX - canvas.offsetLeft, offsetY: touch.pageY - canvas.offsetTop };
        } else {
            return { offsetX: event.offsetX, offsetY: event.offsetY };
        }
    }

    

    function getChargeFromEvent(event) {
        const { offsetX, offsetY } = coordinatesFromMouseOrTouch(event);
        let allowRadius = null;
        if (event.touches) {
            allowRadius = event.touches[0].radiusX + chargeSize;
        } 
        return getChargeFromPoint(charges, offsetX, offsetY, allowRadius);
    }

    function mouseOrTouchDown(event) {
        deselectCharge();
        draggingCharge = getChargeFromEvent(event);
        if (draggingCharge) {
            event.preventDefault();
            const { offsetX, offsetY } = coordinatesFromMouseOrTouch(event);
            draggingOffsetX = offsetX - draggingCharge.x;
            draggingOffsetY = offsetY - draggingCharge.y;
            originalChargeX = draggingCharge.x;
            originalChargeY = draggingCharge.y;
        }
        return false;
    }

    function mouseOrTouchUp(event) {
        if (!draggingCharge) return;
        event.preventDefault();
        if (!(Math.abs(originalChargeX - draggingCharge.x) > 3 || Math.abs(originalChargeY - draggingCharge.y) > 3)) {
            selectCharge(draggingCharge);
        }
        draggingCharge = null;
    }

    function mouseOrTouchMove(event) {
        const { offsetX, offsetY } = coordinatesFromMouseOrTouch(event);
        if (draggingCharge) {
            event.preventDefault();
            draggingCharge.x = offsetX - draggingOffsetX;
            draggingCharge.y = offsetY - draggingOffsetY;
            if(!dynamic)    
                window.requestAnimationFrame(drawVectorField);
        }
    }

    function mouseLeaveOrTouchCancel(event) {
        draggingCharge = null;
    }

    canvas.addEventListener('mousedown', mouseOrTouchDown);
    canvas.addEventListener('mouseup', mouseOrTouchUp);
    canvas.addEventListener('mousemove', mouseOrTouchMove);
    canvas.addEventListener('mouseleave', mouseLeaveOrTouchCancel);

    canvas.addEventListener('touchstart', mouseOrTouchDown);
    canvas.addEventListener('touchend', mouseOrTouchUp);
    canvas.addEventListener('touchmove', mouseOrTouchMove);
    canvas.addEventListener('touchcancel', mouseLeaveOrTouchCancel);

    canvas.addEventListener('dblclick', (event) => {
        let existingCharge = getChargeFromEvent(event);
        if (existingCharge!==null) {
            deselectCharge();
            deleteCharge(existingCharge);
            drawVectorField();
        } else {
            const { offsetX, offsetY } = event;
            deselectCharge();
            addCharge(offsetX, offsetY, nextChargeStrength);
            selectCharge(charges[charges.length - 1]);
        }
        return false;
    });

    addPositiveChargeBtn.addEventListener('click', () => {
        addCharge(canvas.width / (2*dpr), canvas.height / (2*dpr), 1);
    });

    addNegativeChargeBtn.addEventListener('click', () => {
        addCharge(canvas.width / (2*dpr), canvas.height / (2*dpr), -1);
    });

    /*
    fieldStrength.addEventListener('input', () => {
        const strength = parseFloat(fieldStrength.value);
        field.set_ey_field_strength(strength);
        drawVectorField();
    });
    */

    drawVectorField();
}

window.addEventListener('load', main);