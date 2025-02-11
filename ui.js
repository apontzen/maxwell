import init, { compute_field_electrostatic_direct_to_buffer, compute_field_magnetostatic_direct_to_buffer, 
    compute_field_electrostatic_per_charge_direct_to_buffer, compute_field_magnetostatic_per_charge_direct_to_buffer, 
    compute_electric_field_dynamic_to_buffer, init_panic_hook, FieldConfiguration, 
    compute_forces_electrostatic, compute_forces_magnetostatic} from './maxwell/out/maxwell.js';

import { draw, getChargeFromPoint } from './draw.js';

const DYNAMIC_RUNTIME_SECONDS = 10;
const DYNAMIC_RUNTIME_MS = DYNAMIC_RUNTIME_SECONDS * 1000;

let isInitialized = false;
let isInitializing = false;

async function initialize() {
    if(isInitialized) return;
    if(isInitializing) {
        await new Promise(resolve => {
            const interval = setInterval(() => {
                if (isInitialized) {
                    clearInterval(interval);
                    resolve();
                }
            }, 10);
        });
    } else {
        isInitializing = true;
        await init();
        init_panic_hook();
        isInitialized = true;
    }
}
export async function main(params) {

    const {canvas, addPositiveChargeButton, clearChargesButton, solverDropdown,
        potentialControlsDiv, potentialCheckbox, copyJsonButton, pasteJsonButton, 
        chargeOrCurrentSpans, chargePropertiesDiv, startingState, fieldlinesControlsDiv,
        fieldlinesCheckbox, allowEditChargeStrength, allowAddDeleteCharge, 
        perChargeControlDiv, perChargeCheckbox, forcesCheckbox, forcesControlsDiv,
        dipoleControlsDiv, dipoleCheckbox, allowGaussianSurface} = params;
        
    await initialize();

    const ctx = canvas.getContext('2d');
    let rect = canvas.getBoundingClientRect();
    
    let charges = [];
    let draggingCharge = null;
    let draggingDipoleCom = false;
    let draggingOffsetX = 0;
    let draggingOffsetY = 0;

    let drawingGaussianSurface = false;
    let gaussianSurfacePoints = [];
    
    let uniformElecX = 0.0;
    let uniformElecY = 0.0;

    let nextChargeStrength = 1;
    let cic_max_resolution = 128;
    let cic_target_logical_pixel_per_cell = 6;
    let cic_resolution;

    const dpr = window.devicePixelRatio || 1;

    let computeField = compute_field_electrostatic_direct_to_buffer;
    let computeForces = compute_forces_electrostatic;

    let field = null;

    let dynamic = false;

    let showPotential = false;
    let perCharge = false;
    let solver = 'electrostatic_direct';
    let dipoleMode = false;
    let extraDrawInfo = {};
    
    let plotType = 'quiver';

    let animation_request_id = null;
    let last_time = null, last_interaction_time = null;


    function updateSize() {

        // make the DOM element smaller if it's exceeding the width of the device
        canvas.style.width = '90%';
        canvas.style.height = '90%';
        canvas.style.maxWidth = '800px';
        canvas.style.maxHeight = '600px';
        
        // Store the original CSS dimensions.
        rect = canvas.getBoundingClientRect();

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
    }

    
    function animationIsRunning() {
        return animation_request_id !== null && dynamic;
    }

    function continueAnimation() {
        if(dynamic && last_interaction_time !== null && (performance.now() - last_interaction_time) < DYNAMIC_RUNTIME_MS) {
            animation_request_id = window.requestAnimationFrame(tickField);
        }
    }

    function postUserInteraction() {
        last_interaction_time = performance.now();
        saveState();
        if(!animationIsRunning())
            animation_request_id = window.requestAnimationFrame(tickField);
    }

    function solverSupport(solver) {
        const allowPotential = solver === 'electrostatic_direct';
        const allowForces = solver === 'electrostatic_direct' || solver === 'magnetostatic_direct';
        const allowPerCharge = allowForces;
        return {allowPotential, allowForces, allowPerCharge};
    }

    function uiFromState() {
        const {allowPotential, allowForces, allowPerCharge} = solverSupport(solver);

        if(forcesControlsDiv !== null) {
            if(allowForces) {
                forcesControlsDiv.style.display = 'inline';
            } else {
                forcesControlsDiv.style.display = 'none';
            }
        }

        if(potentialControlsDiv !== null) {
            
            if(allowPotential) {
                potentialControlsDiv.style.display = 'inline';

                if(fieldlinesCheckbox !== null)
                    fieldlinesCheckbox.checked = (plotType === 'fieldline');

            } else {
                potentialControlsDiv.style.display = 'none';
            }
        }

        if(perChargeControlDiv !== null) {
            if(allowPerCharge) {
                perChargeControlDiv.style.display = 'inline';
                if(perChargeCheckbox !== null)
                    perChargeCheckbox.checked = perCharge;
            } else {
                perChargeControlDiv.style.display = 'none';
            }
        }

        if(fieldlinesControlsDiv !== null) {
            const allowFieldlines = solver === 'electrostatic_direct' || solver === 'magnetostatic_direct';

            if(allowFieldlines) {
                fieldlinesControlsDiv.style.display = 'inline';
                if(fieldlinesCheckbox !== null)
                    fieldlinesCheckbox.checked = (plotType === 'fieldline');
            } else {
                fieldlinesControlsDiv.style.display = 'none';
                plotType = 'quiver';
            }
        }

        if(solverDropdown !== null)
            solverDropdown.value = solver;

        if(potentialCheckbox !== null)
            potentialCheckbox.checked = showPotential;

        if(dipoleCheckbox !== null) {
            dipoleCheckbox.checked = dipoleMode;
        }

        

        if(solver === 'magnetostatic_direct') {
            updateChargeOrCurrentLabel('Current');
        } else {
            updateChargeOrCurrentLabel('Charge');
        }

    }

    function stateFromUi() {
        
        if(animation_request_id!==null)
            window.cancelAnimationFrame(animation_request_id);

        animation_request_id = null;
        last_time = null;

        if(solverDropdown !== null) 
            solver = solverDropdown.value;

        const solverSupportInfo = solverSupport(solver);

        if(potentialCheckbox !== null)
            showPotential = potentialCheckbox.checked;

        if(fieldlinesCheckbox !== null)
            plotType = fieldlinesCheckbox.checked ? 'fieldline' : 'quiver';

        if(perChargeCheckbox !== null)
            perCharge = solverSupportInfo.allowPerCharge && perChargeCheckbox.checked;

        if(dipoleCheckbox !== null) {
            dipoleMode = dipoleCheckbox.checked;
        }
    
        updateSolverType(true);

    }



    function updateSolverType(is_user_update = false) {
        cic_resolution = Math.min(cic_max_resolution, Math.max(1, Math.floor(rect.width / cic_target_logical_pixel_per_cell)));

        field = new FieldConfiguration(rect.width, rect.height, cic_resolution, cic_resolution);
        
        dynamic = solver === 'dynamic';

        computeForces = null;

        if (solver === 'electrostatic_direct') {
            if(perCharge) {
                computeField = compute_field_electrostatic_per_charge_direct_to_buffer;
            } else {
                computeField = compute_field_electrostatic_direct_to_buffer;
            }
            computeForces = compute_forces_electrostatic;
        } else if (solver === 'electrostatic_fourier' || solver === 'dynamic') {
            computeField = compute_electric_field_dynamic_to_buffer;
        } else if (solver === 'magnetostatic_direct') { 
            if(perCharge) {
                computeField = compute_field_magnetostatic_per_charge_direct_to_buffer;
            } else { 
                computeField = compute_field_magnetostatic_direct_to_buffer;
            }
            computeForces = compute_forces_magnetostatic;
        } else {
            console.error('Unknown solver type');
        }

        uiFromState(); // e.g. updates visibility of potential checkbox

        if(is_user_update) {
            postUserInteraction();
        } else {
            drawVectorField();
        }
    }

    if(fieldlinesCheckbox!==null)
        fieldlinesCheckbox.addEventListener('change', stateFromUi);
    if(forcesCheckbox!==null)
        forcesCheckbox.addEventListener('change', stateFromUi);
    if(solverDropdown !== null)
        solverDropdown.addEventListener('change', stateFromUi);
    if(potentialCheckbox !== null)
        potentialCheckbox.addEventListener('change', stateFromUi);
    if(perChargeCheckbox !== null)
        perChargeCheckbox.addEventListener('change', stateFromUi);
    if(dipoleCheckbox !== null)
        dipoleCheckbox.addEventListener('change', stateFromUi);

    updateSize();
    updateSolverType();
    drawVectorField();
    
    

    function stateToJson() {
        let charges_normalized_coordinates = charges.map(charge => {
            const normalizedX = charge.x / rect.width;
            const normalizedY = charge.y / rect.height;
            return { ...charge, x: normalizedX, y: normalizedY };
        });

        return JSON.stringify({
            charges: charges_normalized_coordinates,
            solver: solver,
            showPotential: showPotential,
            plotType: plotType,
            perCharge: perCharge,
            // for information about the scaling of the uniform field, see jsonToState
            uniformElecX: uniformElecX * rect.height, 
            uniformElecY: uniformElecY * rect.width,
            dipoleMode: dipoleMode,
            extraDrawInfo: extraDrawInfo
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
            
            showPotential = state.showPotential;

            if(state.plotType)
                plotType = state.plotType;

            // The below scaling looks weird but it is motivated by allowing a line
            // of a fixed number of charges to cancel out the uniform field. This
            // means the electric field needs to scale inversely with the perpendicular
            // dimension.
            if(state.uniformElecX)
                uniformElecX = state.uniformElecX / rect.height;
            else
                uniformElecX = 0.0;

            if(state.uniformElecY)
                uniformElecY = state.uniformElecY / rect.width;
            else
                uniformElecY = 0.0;

            solver = state.solver;

            perCharge = state.perCharge;

            dipoleMode = state.dipoleMode;

            if (state.extraDrawInfo !== undefined)
                extraDrawInfo = state.extraDrawInfo;

            updateSolverType();

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
            console.log("Restoring saved state", savedState);
            jsonToState(savedState);
        } else {
            console.log("No saved state found");
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

        const dx = rect.width / cic_resolution;

        // Try to get forward in time by dt, but abandon if it takes too long (more than 20ms)
        const performance_time_start = performance.now();
        const max_time_ms = 20.0;

        while (dt > 0.0 && performance.now() - performance_time_start < max_time_ms) {
            const step = Math.min(0.1 * dx , dt);
            field.tick(step);
            dt -= step;
        }
 
        last_time = time_now;
        drawVectorField(false);
 

        // only request the next frame if we're still solving dynamically
        

    }

    function closeLoop(points) {
        if(points.length < 2) return [];
        return [...points, points[0]];
    }

    function drawVectorField() {
        let charges_no_test_charges = charges.filter(c => !c.isTestCharge);
        field.set_charges(charges_no_test_charges);
        field.set_uniform_field(uniformElecX, uniformElecY);

        if(!dynamic) {
            // if dynamic, the fields are being updated in real time and we don't want to reset them
            // if not dynamic, we need to reset the fields to trigger a recompute
            field.reset_fields();
        }


        let forces = null;

        if(forcesCheckbox!==null && forcesCheckbox.checked && computeForces !== null) {
            forces = computeForces(field);
        }

        extraDrawInfo['gaussianSurfacePoints'] = closeLoop(gaussianSurfacePoints);

        draw(ctx, rect, charges, field, plotType, computeField, 
            computeField === compute_field_electrostatic_direct_to_buffer && showPotential, null,
            forces, dipoleMode, extraDrawInfo);

        continueAnimation();

    }

    function updateChargeOrCurrentLabel(charge_or_current) {
        // find all spans with class charge_or_current and update their contents to the provided string
        if(chargeOrCurrentSpans!==null)
            chargeOrCurrentSpans.forEach(span => {
                span.textContent = charge_or_current;
            });
    }

    
    
    let selectedCharge = null;

    let chargeInput, chargeValue, deleteChargeBtn, testChargeTickbox;

    if(chargePropertiesDiv!==null) {
        chargeInput = chargePropertiesDiv.querySelector('.charge');
        chargeValue = chargePropertiesDiv.querySelector('.chargeValue');
        deleteChargeBtn = chargePropertiesDiv.querySelector('.deleteCharge');
        testChargeTickbox = chargePropertiesDiv.querySelector('.testCharge');

        chargeInput.addEventListener('input', () => {
            selectedCharge.charge = parseInt(chargeInput.value);
            chargeValue.textContent = chargeInput.value;
            drawVectorField();
        });

        testChargeTickbox.addEventListener('change', () => {
            selectedCharge.isTestCharge = testChargeTickbox.checked;
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
        testChargeTickbox.checked = charge.isTestCharge;
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
            const canvasRect = canvas.getBoundingClientRect();
            const offsetTop  = canvasRect.top + window.scrollY;
            const offsetLeft = canvasRect.left + window.scrollX;

            return { offsetX: touch.pageX - offsetLeft, offsetY: touch.pageY - offsetTop };
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
        return getChargeFromPoint(charges, offsetX, offsetY, allowRadius, true, null, dipoleMode);
    }

    function mouseOrTouchDown(event) {
        deselectCharge();
        draggingCharge = getChargeFromEvent(event);
        if (Array.isArray(draggingCharge)) {
            draggingDipoleCom = true;
            draggingCharge = draggingCharge[0];
        } else {
            draggingDipoleCom = false;
        }
        
        if (draggingCharge) {
            event.preventDefault();
            const { offsetX, offsetY } = coordinatesFromMouseOrTouch(event);
            draggingOffsetX = offsetX - draggingCharge.x;
            draggingOffsetY = offsetY - draggingCharge.y;
            originalChargeX = draggingCharge.x;
            originalChargeY = draggingCharge.y;
        } else if(allowGaussianSurface) {
            drawingGaussianSurface = true;
            event.preventDefault();
            const { offsetX, offsetY } = coordinatesFromMouseOrTouch(event);
            gaussianSurfacePoints = [{x: offsetX, y: offsetY}];
        }
        return false;
    }

    function mouseOrTouchUp(event) {
        if(drawingGaussianSurface) {
            drawingGaussianSurface = false;
            event.preventDefault();
        } else if(draggingCharge) {
            event.preventDefault();
            if (!(Math.abs(originalChargeX - draggingCharge.x) > 3 || Math.abs(originalChargeY - draggingCharge.y) > 3)) {
                selectCharge(draggingCharge);
            }
            draggingCharge = null;
        }
    }

    function mouseOrTouchMove(event) {
        const { offsetX, offsetY } = coordinatesFromMouseOrTouch(event);
        if(drawingGaussianSurface) {
            const distanceFromLast = Math.sqrt((offsetX - gaussianSurfacePoints[gaussianSurfacePoints.length - 1].x) ** 2 + (offsetY - gaussianSurfacePoints[gaussianSurfacePoints.length - 1].y) ** 2);
            if(distanceFromLast > 10) {
                gaussianSurfacePoints.push({x: offsetX, y: offsetY});
                postUserInteraction();
            }
        } else if (draggingCharge) {
            event.preventDefault();
            let otherCharge, deltaX, deltaY, comX, comY;

           
            if(dipoleMode && draggingCharge.dipoleWith !== undefined) {
                otherCharge = charges.find(c => c.id === draggingCharge.dipoleWith);
                deltaX = otherCharge.x - draggingCharge.x;
                deltaY = otherCharge.y - draggingCharge.y;
                comX = (otherCharge.x + draggingCharge.x) / 2;
                comY = (otherCharge.y + draggingCharge.y) / 2;
            }

            draggingCharge.x = offsetX - draggingOffsetX;
            draggingCharge.y = offsetY - draggingOffsetY;


            if(otherCharge !== undefined) {
                if(!draggingDipoleCom) {
                    let newDeltaX = draggingCharge.x - comX;
                    let newDeltaY = draggingCharge.y - comY;

                    // normalize to same distance as original delta:
                    let scale = Math.sqrt(deltaX * deltaX + deltaY * deltaY) / Math.sqrt(newDeltaX * newDeltaX + newDeltaY * newDeltaY);
                    newDeltaX *= scale;
                    newDeltaY *= scale;

                    draggingCharge.x = comX + newDeltaX/2;
                    draggingCharge.y = comY + newDeltaY/2;
                    otherCharge.x = comX - newDeltaX/2;
                    otherCharge.y = comY - newDeltaY/2;
                } else {
                    otherCharge.x = draggingCharge.x + deltaX;
                    otherCharge.y = draggingCharge.y + deltaY;
                }
            }
            
            postUserInteraction();
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

    window.addEventListener('resize', () => {
        const state = stateToJson();
        updateSize();
        jsonToState(state);
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

}

export function initialize_on_existing_dom() {
    window.addEventListener('load', ()=> {
        const canvas = document.getElementById('vectorFieldCanvas');
        const addPositiveChargeButton = document.getElementById('addPositiveCharge');
        const clearChargesButton = document.getElementById('clearCharges');
        const solverDropdown = document.getElementById('solver');
        const potentialControlsDiv = document.getElementById('show-potential-control');
        const fieldlinesControlsDiv = document.getElementById('fieldlines-control');
        const potentialCheckbox = document.getElementById('potential');
        const fieldlinesCheckbox = document.getElementById('fieldlines');
        const copyJsonButton = document.getElementById('copyJson');
        const pasteJsonButton = document.getElementById('pasteJson');
        const chargeOrCurrentSpans = document.querySelectorAll('.charge-or-current');
        const chargePropertiesDiv = document.querySelector('.charge-properties');
        const perChargeControlDiv = document.getElementById('per-charge-control');
        const perChargeCheckbox = document.getElementById('per-charge');
        const forcesCheckbox = document.getElementById('forces');
        const forcesControlsDiv = document.getElementById('forces-control');
        const dipoleControlsDiv = document.getElementById('dipole-control');
        const dipoleCheckbox = document.getElementById('dipole');
        const allowGaussianSurface = false;
        const startingState = null;

        const allowEditChargeStrength = true;
        const allowAddDeleteCharge = true;

        const params = {
            canvas, addPositiveChargeButton, clearChargesButton, solverDropdown, potentialControlsDiv,
            fieldlinesControlsDiv, fieldlinesCheckbox, 
            potentialCheckbox, copyJsonButton, pasteJsonButton, chargeOrCurrentSpans, chargePropertiesDiv,
            startingState, allowEditChargeStrength, allowAddDeleteCharge, perChargeControlDiv, perChargeCheckbox,
            forcesCheckbox, forcesControlsDiv, dipoleControlsDiv, dipoleCheckbox, allowGaussianSurface
        }
        main(params);

    });
}

function addCheckbox(className, labelText, parentElement) {
    parentElement.style.display = 'inline-block';
    
    const controlDiv = document.createElement('div');
    controlDiv.className = className;

    const checkbox = document.createElement('input');
    checkbox.type = 'checkbox';
    checkbox.checked = false;
    checkbox.id = className + '-' + Math.random().toString(36).substring(7);
    controlDiv.appendChild(checkbox);

    const label = document.createElement('label');
    label.htmlFor = checkbox.id;
    label.textContent = labelText;
    controlDiv.appendChild(label);

    parentElement.appendChild(controlDiv);

    return {div: controlDiv, checkbox};
}

export function embed() {
    const meme = document.createElement('meme-embed');
    window.addEventListener('load', () => {
        const all_memes = document.querySelectorAll('meme-embed');
        all_memes.forEach(meme =>  {
            meme.style.position ='relative';
            const startingState = meme.getAttribute('meme');

            let fieldlinesControlsDiv = null;
            let fieldlinesCheckbox = null;
            let perChargeControlDiv = null;
            let perChargeCheckbox = null;
            let forcesControlsDiv = null;
            let forcesCheckbox = null;
            let dipoleControlsDiv = null;
            let dipoleCheckbox = null;

            if(meme.getAttribute('fieldlines-checkbox')) {
                ({div: fieldlinesControlsDiv, checkbox: fieldlinesCheckbox} = addCheckbox('fieldlines-controls', 'Show field lines', meme));
            }

            if(meme.getAttribute('percharge-checkbox')) {
                let label = 'Per charge';
                if (startingState && JSON.parse(startingState).solver === 'magnetostatic_direct') {
                    label = 'Per current';
                }
                ({div: perChargeControlDiv, checkbox: perChargeCheckbox} = addCheckbox('percharge-controls', label, meme));
            }

            if(meme.getAttribute('forces-checkbox')) {
                ({div: forcesControlsDiv, checkbox: forcesCheckbox} = addCheckbox('forces-controls', 'Show forces', meme));
            }

            if(meme.getAttribute('dipole-checkbox')) {
                ({div: dipoleControlsDiv, checkbox: dipoleCheckbox} = addCheckbox('dipole-controls', 'Rigid dipole', meme));
            }



            const canvas = meme.appendChild(document.createElement('canvas'));


            const params = {
                canvas, addPositiveChargeButton: null, clearChargesButton: null, solverDropdown: null, potentialControlsDiv: null,
                potentialCheckbox: null, copyJsonButton: null, pasteJsonButton: null, chargeOrCurrentSpans: null, chargePropertiesDiv: null,
                startingState: startingState, allowEditChargeStrength: false, allowAddDeleteCharge: false,
                fieldlinesControlsDiv, fieldlinesCheckbox, perChargeControlDiv, perChargeCheckbox, forcesCheckbox, forcesControlsDiv,
                dipoleControlsDiv, dipoleCheckbox,
                allowGaussianSurface: meme.getAttribute('gaussian-surface') !== null
            }
            main(params);
        });
    });
}