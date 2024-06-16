import init, { compute_field_electrostatic_direct, compute_field_electrostatic_fourier, init_panic_hook, FieldConfiguration } from './maxwell/out/maxwell.js';

export async function main() {
    await init();

    init_panic_hook();


    const canvas = document.getElementById('vectorFieldCanvas');
    const ctx = canvas.getContext('2d');
    const plotTypeSelect = document.getElementById('plotType');
    const addPositiveChargeBtn = document.getElementById('addPositiveCharge');
    const addNegativeChargeBtn = document.getElementById('addNegativeCharge');

    const chargeSize = 10;
    const maxArrowLength = 40;

    let charges = [];
    let draggingCharge = null;
    let draggingOffsetX = 0;
    let draggingOffsetY = 0;

    let cic_resolution = 20;

    const dpr = window.devicePixelRatio || 1;

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
    
    function updateSolverType() {
        // get value of solver
        const solverType = document.getElementById('solver').value;
        if (solverType === 'electrostatic_direct') {
            computeField = compute_field_electrostatic_direct;
        } else if (solverType === 'electrostatic_fourier') {
            computeField = compute_field_electrostatic_fourier;
        } else {
            console.error('Unknown solver type');
        }
    }

    updateSolverType();

    document.getElementById('solver').addEventListener('change', () => {
        updateSolverType();
        drawVectorField();
    });

    function drawVectorField() {
        const plotType = plotTypeSelect.value;
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        const vectors = generateVectors();
        
        if (plotType === 'quiver') {
            drawQuiverPlot(vectors);
        } else if (plotType === 'streamline') {
            drawStreamlinePlot();
        }

        drawCharges();
    }

    function generateVectors() {
        const vectors = [];
        const step = 20;
        const field = new FieldConfiguration(canvas.width, canvas.height, cic_resolution, cic_resolution);
        field.set_charges(charges);

        for (let x = 0; x < canvas.width; x += step) {
            for (let y = 0; y < canvas.height; y += step) {
                // exclude the vector if it's within step distance from any charge
                if (charges.some(charge => {
                    const dx = x - charge.x;
                    const dy = y - charge.y;
                    return Math.sqrt(dx * dx + dy * dy) < step;
                })) continue;
                const vector = computeField(field, x, y);
                vectors.push({x, y, u: vector[0], v: vector[1]});
            }
        }
        
        return vectors;
    }

    function drawQuiverPlot(vectors) {
        vectors.forEach(({x, y, u, v}) => {
            drawArrow(x, y, u, v);
        });
    }

    function drawStreamlinePlot() {
        const numStreamlines = 30;
        const maxLength = 1000;
        const stepSize = 2;
        for (let i = 0; i < numStreamlines; i++) {
            let x = Math.random() * canvas.width;
            let y = Math.random() * canvas.height;
            ctx.beginPath();
            ctx.moveTo(x, y);
            for (let j = 0; j < maxLength; j++) {
                const {u, v} = computeField(x, y);
                const speed = Math.sqrt(u * u + v * v);
                if (speed < 0.1) break;
                x += (u / speed) * stepSize;
                y += (v / speed) * stepSize;
                if (x < 0 || x >= canvas.width || y < 0 || y >= canvas.height) break;
                ctx.lineTo(x, y);
            }
            ctx.stroke();
        }
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
    let selectedCharge = null;

    chargeInput.addEventListener('input', () => {
        selectedCharge.charge = parseInt(chargeInput.value);
        drawVectorField();
    });

    const deleteChargeBtn = document.getElementById('deleteCharge');
    deleteChargeBtn.addEventListener('click', () => {
        charges = charges.filter(c => c !== selectedCharge);
        deselectCharge();
        drawVectorField();
    });


    function selectCharge(charge) {

        chargeInput.value = charge.charge;
        selectedCharge = charge;

        const chargeProperties = document.querySelector('.charge-properties');
        chargeProperties.style.display = 'block';
    }

    function deselectCharge() {
        const chargeProperties = document.querySelector('.charge-properties');
        chargeProperties.style.display = 'none';
        selectedCharge = null;
    }

    function drawCharges() {
        charges.forEach(charge => {
            ctx.beginPath();
            ctx.arc(charge.x, charge.y, chargeSize, 0, 2 * Math.PI, false);
            ctx.fillStyle = charge.charge > 0 ? 'red' : 'blue';
            ctx.fill();
            ctx.stroke();
        });
    }

    function addCharge(x, y, charge) {
        charges.push({x, y, charge});
        drawVectorField();
    }

    canvas.addEventListener('mousedown', (event) => {
        const {offsetX, offsetY} = event;
        for (const charge of charges) {
            const dx = offsetX - charge.x;
            const dy = offsetY - charge.y;
            if (Math.sqrt(dx * dx + dy * dy) < 10) {
                draggingCharge = charge;
                draggingOffsetX = dx;
                draggingOffsetY = dy;
                break;
            }
        }
    });

    canvas.addEventListener('mousemove', (event) => {
        if (draggingCharge) {
            draggingCharge.x = event.offsetX - draggingOffsetX;
            draggingCharge.y = event.offsetY - draggingOffsetY;
            drawVectorField();
        }
    });

    canvas.addEventListener('mouseup', () => {
        draggingCharge = null;
    });

    canvas.addEventListener('click', (event) => {
        const {offsetX, offsetY} = event;
        const charge = charges.find(charge => {
            const dx = offsetX - charge.x;
            const dy = offsetY - charge.y;
            return Math.sqrt(dx * dx + dy * dy) < 10;
        });
        deselectCharge();
        if (charge) {
            selectCharge(charge);
        } 
    });

    canvas.addEventListener('mouseleave', () => {
        draggingCharge = null;
    });

    addPositiveChargeBtn.addEventListener('click', () => {
        addCharge(canvas.width / (2*dpr), canvas.height / (2*dpr), 1);
    });

    addNegativeChargeBtn.addEventListener('click', () => {
        addCharge(canvas.width / (2*dpr), canvas.height / (2*dpr), -1);
    });

    plotTypeSelect.addEventListener('change', drawVectorField);

    drawVectorField();
}

window.addEventListener('load', main);