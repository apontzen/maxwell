import init, { compute_field_electrostatic_direct_to_buffer, compute_field_magnetostatic_direct_to_buffer,
    compute_electric_field_dynamic_to_buffer, init_panic_hook, FieldConfiguration } from './maxwell/out/maxwell.js';

import { draw, getChargeFromPoint } from './draw.js';

let isInitialized = false;




export async function main(params) {

    const {canvas, addPositiveChargeButton, clearChargesButton, solverDropdown,
        potentialControlsDiv, potentialCheckbox, copyJsonButton, pasteJsonButton, 
        chargeOrCurrentSpans, chargePropertiesDiv, startingState,
        allowEditChargeStrength, allowAddDeleteCharge} = params;
        
    if(!isInitialized) {
        await init();
        init_panic_hook();
        isInitialized = true;
    }

    const ctx = canvas.getContext('2d');
    

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

    let computeField = compute_field_electrostatic_direct_to_buffer;
    let field = null;

    let dynamic = false;

    let showPotential = false;
    let solver = null;
    
    let plotType = 'quiver';

    function updateSolverType(solverType = null) {
        field = new FieldConfiguration(rect.width, rect.height, cic_resolution, cic_resolution);

        if (solverDropdown === null && solverType === null)
            return;

        if (solverType === null) {
            solverType = solverDropdown.value;
        }

        solver = solverType;
        
        let allowPotential = false;

        dynamic = solverType === 'dynamic';

        plotType = 'quiver';

        if (solverType.endsWith('_fieldline')) {
            plotType = 'fieldline';
        }

        if (solverType === 'electrostatic_direct' || solverType === 'electrostatic_direct_fieldline') {
            updateChargeOrCurrentLabel('Charge');
            computeField = compute_field_electrostatic_direct_to_buffer;
            allowPotential = true;
        } else if (solverType === 'electrostatic_fourier' || solverType === 'dynamic') {
            updateChargeOrCurrentLabel('Charge');
            computeField = compute_electric_field_dynamic_to_buffer;
        } else if (solverType === 'magnetostatic_direct' || solverType === 'magnetostatic_direct_fieldline') { 
            updateChargeOrCurrentLabel('Current');
            computeField = compute_field_magnetostatic_direct_to_buffer;
        } else {
            console.error('Unknown solver type');
        }

        if(potentialControlsDiv !== null) {
            if(allowPotential) {
                potentialControlsDiv.style.display = 'inline';
            } else {
                potentialControlsDiv.style.display = 'none';
            }
        }

        if(dynamic)
            animation_request_id = window.requestAnimationFrame(tickField);
    }

    updateSolverType();

    let animation_request_id = null;
    let last_time = null;

    if(solverDropdown !== null)
        solverDropdown.addEventListener('change', () => {
            if(animation_request_id!==null)
                window.cancelAnimationFrame(animation_request_id);
            animation_request_id = null;
            last_time = null;
            updateSolverType();
            drawVectorField();
        });

    if(potentialCheckbox !== null)
        potentialCheckbox.addEventListener('change', () => {
            showPotential = potentialCheckbox.checked;
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
            solver: solver,
            show_potential: showPotential
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

            if(solverDropdown !== null)
                solverDropdown.value = state.solver;

            if(potentialCheckbox !== null)
                potentialCheckbox.checked = state.show_potential;
            
            showPotential = state.show_potential;
            updateSolverType(state.solver);

            return true;
        } else {
            return false;
        }
    }

    if(copyJsonButton !== null)
        copyJsonButton.addEventListener('click', () => {
            const json = stateToJson();
            navigator.clipboard.writeText(json);
        });

    if(pasteJsonButton !== null)
        pasteJsonButton.addEventListener('click', async () => {
            const json = await navigator.clipboard.readText();
            if(jsonToState(json)) {
                drawVectorField();
            }
        });

    if(startingState) {
        jsonToState(startingState);
    } else { 
        // Load the state from local storage if it exists
        const savedState = localStorage.getItem('maxwell_state');
        if (savedState) {
            jsonToState(savedState);
        }
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
        field.set_charges(charges);

        if(!dynamic) {
            field.reset_fields();
            // if dynamic, the fields are being updated in real time and we don't want to reset them
        }

        if(user_update) {
            saveState();
        }

        draw(ctx, rect, charges, field, plotType, computeField, 
            computeField === compute_field_electrostatic_direct_to_buffer && showPotential);



    }

    function updateChargeOrCurrentLabel(charge_or_current) {
        // find all spans with class charge_or_current and update their contents to the provided string
        if(chargeOrCurrentSpans!==null)
            chargeOrCurrentSpans.forEach(span => {
                span.textContent = charge_or_current;
            });
    }

    
    
    let selectedCharge = null;

    let chargeInput, chargeValue, deleteChargeBtn;

    if(chargePropertiesDiv!==null) {
        chargeInput = chargePropertiesDiv.querySelector('.charge');
        chargeValue = chargePropertiesDiv.querySelector('.chargeValue');
        deleteChargeBtn = chargePropertiesDiv.querySelector('.deleteCharge');

        chargeInput.addEventListener('input', () => {
            selectedCharge.charge = parseInt(chargeInput.value);
            chargeValue.textContent = chargeInput.value;
            drawVectorField();
        });

        function deleteCharge(charge) {
            if(!allowAddDeleteCharge) return;
            charges = charges.filter(c => c !== charge);
            field.set_charges(charges);
        }

        deleteChargeBtn.addEventListener('click', () => {
            deleteCharge(selectedCharge);
            deselectCharge();
            drawVectorField();
        });
    }

    function selectCharge(charge) {

        if(!allowEditChargeStrength) return;
        if(chargePropertiesDiv===null) return;

        chargeInput.value = charge.charge;
        chargeValue.textContent = chargeInput.value;
        selectedCharge = charge;

        drawVectorField();

        
        chargePropertiesDiv.style.display = 'block';

        chargePropertiesDiv.style.display = 'block';
        chargePropertiesDiv.style.position = 'absolute';
        const canvasRect = canvas.getBoundingClientRect();
        const pageOffsetTop = window.scrollY || document.documentElement.scrollTop || document.body.scrollTop || 0;
        const pageOffsetLeft = window.scrollX || document.documentElement.scrollLeft || document.body.scrollLeft || 0;
        const canvasRectTop = canvasRect.top + pageOffsetTop;
        const canvasRectLeft = canvasRect.left + pageOffsetLeft;
        chargePropertiesDiv.style.top = `${canvasRectTop + charge.y - chargePropertiesDiv.offsetHeight / 2}px`;

        // Check if charge is on the right side of the screen
        if (charge.x > rect.width / 2) {
            chargePropertiesDiv.style.left = `${canvasRectLeft + charge.x - chargePropertiesDiv.offsetWidth - 20}px`;
            chargePropertiesDiv.classList.add('point-right');
            chargePropertiesDiv.classList.remove('point-left');
        } else {
            chargePropertiesDiv.style.left = `${canvasRectLeft + charge.x + 20}px`;
            chargePropertiesDiv.classList.add('point-left');
            chargePropertiesDiv.classList.remove('point-right');
        }
        // if the chargeProperties div is off the screen, move it back on
        if (chargePropertiesDiv.offsetLeft < 0) {
            chargePropertiesDiv.style.left = '0';
        }
        if (chargePropertiesDiv.offsetLeft + chargePropertiesDiv.offsetWidth > window.innerWidth) {
            chargePropertiesDiv.style.left = `${window.innerWidth - chargePropertiesDiv.offsetWidth}px`;
        }
        chargePropertiesDiv.style.zIndex = '1';

        
    }

    function deselectCharge() {
    
        if (selectedCharge !== null) {
            nextChargeStrength = selectedCharge.charge;         
        }

        selectedCharge = null;

        if(chargePropertiesDiv!==null) 
            chargePropertiesDiv.style.display = 'none';

        drawVectorField();
    }


    let charge_id = 0;

    function addCharge(x, y, charge) {
        if(!allowAddDeleteCharge) return;
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
        const buffer = new Float64Array(2);
        for (let x = 0; x < rect.width; x += 10) {
            computeField(field, x, y, buffer);
            x_values.push(x);
            u_values.push(buffer[0]);
            v_values.push(buffer[1]);
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
            allowRadius = event.touches[0].radiusX;
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

    if(addPositiveChargeButton!==null) 
        addPositiveChargeButton.addEventListener('click', () => {
            addCharge(canvas.width / (2*dpr), canvas.height / (2*dpr), 1);
        });

    if(clearChargesButton!==null)
        clearChargesButton.addEventListener('click', () => {
            if(!allowAddDeleteCharge) return;
            charges = [];
            deselectCharge();
            drawVectorField();
        });

    drawVectorField();
}

export function initialize_on_existing_dom() {
    window.addEventListener('load', ()=> {
        const canvas = document.getElementById('vectorFieldCanvas');
        const addPositiveChargeButton = document.getElementById('addPositiveCharge');
        const clearChargesButton = document.getElementById('clearCharges');
        const solverDropdown = document.getElementById('solver');
        const potentialControlsDiv = document.getElementById('show-potential-control')
        const potentialCheckbox = document.getElementById('potential');
        const copyJsonButton = document.getElementById('copyJson');
        const pasteJsonButton = document.getElementById('pasteJson');
        const chargeOrCurrentSpans = document.querySelectorAll('.charge-or-current');
        const chargePropertiesDiv = document.querySelector('.charge-properties');

        const startingState = null;

        const allowEditChargeStrength = true;
        const allowAddDeleteCharge = true;

        const params = {
            canvas, addPositiveChargeButton, clearChargesButton, solverDropdown, potentialControlsDiv,
            potentialCheckbox, copyJsonButton, pasteJsonButton, chargeOrCurrentSpans, chargePropertiesDiv,
            startingState, allowEditChargeStrength, allowAddDeleteCharge
        }
        main(params);

    });
}

export function embed() {
    const meme = document.createElement('meme-embed');
    window.addEventListener('load', () => {
        const all_memes = document.querySelectorAll('meme-embed');
        all_memes.forEach(meme => {
            const canvas = meme.appendChild(document.createElement('canvas'));
            const startingState = meme.getAttribute('meme');

            const params = {
                canvas, addPositiveChargeButton: null, clearChargesButton: null, solverDropdown: null, potentialControlsDiv: null,
                potentialCheckbox: null, copyJsonButton: null, pasteJsonButton: null, chargeOrCurrentSpans: null, chargePropertiesDiv: null,
                startingState: startingState, allowEditChargeStrength: false, allowAddDeleteCharge: false
            }
            main(params);
        });
    });
}