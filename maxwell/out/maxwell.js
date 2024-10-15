let wasm;

const cachedTextDecoder = (typeof TextDecoder !== 'undefined' ? new TextDecoder('utf-8', { ignoreBOM: true, fatal: true }) : { decode: () => { throw Error('TextDecoder not available') } } );

if (typeof TextDecoder !== 'undefined') { cachedTextDecoder.decode(); };

let cachedUint8Memory0 = null;

function getUint8Memory0() {
    if (cachedUint8Memory0 === null || cachedUint8Memory0.byteLength === 0) {
        cachedUint8Memory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8Memory0;
}

function getStringFromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return cachedTextDecoder.decode(getUint8Memory0().subarray(ptr, ptr + len));
}

const heap = new Array(128).fill(undefined);

heap.push(undefined, null, true, false);

let heap_next = heap.length;

function addHeapObject(obj) {
    if (heap_next === heap.length) heap.push(heap.length + 1);
    const idx = heap_next;
    heap_next = heap[idx];

    heap[idx] = obj;
    return idx;
}

function getObject(idx) { return heap[idx]; }

function dropObject(idx) {
    if (idx < 132) return;
    heap[idx] = heap_next;
    heap_next = idx;
}

function takeObject(idx) {
    const ret = getObject(idx);
    dropObject(idx);
    return ret;
}

function getArrayU8FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint8Memory0().subarray(ptr / 1, ptr / 1 + len);
}

function isLikeNone(x) {
    return x === undefined || x === null;
}

let cachedFloat64Memory0 = null;

function getFloat64Memory0() {
    if (cachedFloat64Memory0 === null || cachedFloat64Memory0.byteLength === 0) {
        cachedFloat64Memory0 = new Float64Array(wasm.memory.buffer);
    }
    return cachedFloat64Memory0;
}

let cachedInt32Memory0 = null;

function getInt32Memory0() {
    if (cachedInt32Memory0 === null || cachedInt32Memory0.byteLength === 0) {
        cachedInt32Memory0 = new Int32Array(wasm.memory.buffer);
    }
    return cachedInt32Memory0;
}

let WASM_VECTOR_LEN = 0;

const cachedTextEncoder = (typeof TextEncoder !== 'undefined' ? new TextEncoder('utf-8') : { encode: () => { throw Error('TextEncoder not available') } } );

const encodeString = (typeof cachedTextEncoder.encodeInto === 'function'
    ? function (arg, view) {
    return cachedTextEncoder.encodeInto(arg, view);
}
    : function (arg, view) {
    const buf = cachedTextEncoder.encode(arg);
    view.set(buf);
    return {
        read: arg.length,
        written: buf.length
    };
});

function passStringToWasm0(arg, malloc, realloc) {

    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length, 1) >>> 0;
        getUint8Memory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;

    const mem = getUint8Memory0();

    let offset = 0;

    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }

    if (offset !== len) {
        if (offset !== 0) {
            arg = arg.slice(offset);
        }
        ptr = realloc(ptr, len, len = offset + arg.length * 3, 1) >>> 0;
        const view = getUint8Memory0().subarray(ptr + offset, ptr + len);
        const ret = encodeString(arg, view);

        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

function debugString(val) {
    // primitive types
    const type = typeof val;
    if (type == 'number' || type == 'boolean' || val == null) {
        return  `${val}`;
    }
    if (type == 'string') {
        return `"${val}"`;
    }
    if (type == 'symbol') {
        const description = val.description;
        if (description == null) {
            return 'Symbol';
        } else {
            return `Symbol(${description})`;
        }
    }
    if (type == 'function') {
        const name = val.name;
        if (typeof name == 'string' && name.length > 0) {
            return `Function(${name})`;
        } else {
            return 'Function';
        }
    }
    // objects
    if (Array.isArray(val)) {
        const length = val.length;
        let debug = '[';
        if (length > 0) {
            debug += debugString(val[0]);
        }
        for(let i = 1; i < length; i++) {
            debug += ', ' + debugString(val[i]);
        }
        debug += ']';
        return debug;
    }
    // Test for built-in
    const builtInMatches = /\[object ([^\]]+)\]/.exec(toString.call(val));
    let className;
    if (builtInMatches.length > 1) {
        className = builtInMatches[1];
    } else {
        // Failed to match the standard '[object ClassName]'
        return toString.call(val);
    }
    if (className == 'Object') {
        // we're a user defined class or Object
        // JSON.stringify avoids problems with cycles, and is generally much
        // easier than looping through ownProperties of `val`.
        try {
            return 'Object(' + JSON.stringify(val) + ')';
        } catch (_) {
            return 'Object';
        }
    }
    // errors
    if (val instanceof Error) {
        return `${val.name}: ${val.message}\n${val.stack}`;
    }
    // TODO we could test for more things here, like `Set`s and `Map`s.
    return className;
}
/**
*/
export function init_panic_hook() {
    wasm.init_panic_hook();
}

function _assertClass(instance, klass) {
    if (!(instance instanceof klass)) {
        throw new Error(`expected instance of ${klass.name}`);
    }
    return instance.ptr;
}
/**
* @param {FieldConfiguration} field_configuration
* @param {number} x
* @param {number} y
* @returns {number}
*/
export function compute_potential_electrostatic_direct(field_configuration, x, y) {
    _assertClass(field_configuration, FieldConfiguration);
    const ret = wasm.compute_potential_electrostatic_direct(field_configuration.__wbg_ptr, x, y);
    return ret;
}

/**
* @param {FieldConfiguration} field_configuration
* @param {number} x
* @param {number} y
* @returns {Pair}
*/
export function compute_field_electrostatic_direct(field_configuration, x, y) {
    _assertClass(field_configuration, FieldConfiguration);
    const ret = wasm.compute_field_electrostatic_direct(field_configuration.__wbg_ptr, x, y);
    return Pair.__wrap(ret);
}

function passArrayF64ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 8, 8) >>> 0;
    getFloat64Memory0().set(arg, ptr / 8);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}
/**
* @param {FieldConfiguration} field_configuration
* @param {number} x
* @param {number} y
* @param {Float64Array} buffer
*/
export function compute_field_electrostatic_direct_to_buffer(field_configuration, x, y, buffer) {
    _assertClass(field_configuration, FieldConfiguration);
    var ptr0 = passArrayF64ToWasm0(buffer, wasm.__wbindgen_malloc);
    var len0 = WASM_VECTOR_LEN;
    wasm.compute_field_electrostatic_direct_to_buffer(field_configuration.__wbg_ptr, x, y, ptr0, len0, addHeapObject(buffer));
}

/**
* @param {FieldConfiguration} field_configuration
* @param {number} x
* @param {number} y
* @param {Float64Array} buffer
*/
export function compute_field_electrostatic_per_charge_direct_to_buffer(field_configuration, x, y, buffer) {
    _assertClass(field_configuration, FieldConfiguration);
    var ptr0 = passArrayF64ToWasm0(buffer, wasm.__wbindgen_malloc);
    var len0 = WASM_VECTOR_LEN;
    wasm.compute_field_electrostatic_per_charge_direct_to_buffer(field_configuration.__wbg_ptr, x, y, ptr0, len0, addHeapObject(buffer));
}

/**
* @param {FieldConfiguration} field_configuration
* @param {number} x
* @param {number} y
* @returns {Pair}
*/
export function compute_field_magnetostatic_direct(field_configuration, x, y) {
    _assertClass(field_configuration, FieldConfiguration);
    const ret = wasm.compute_field_magnetostatic_direct(field_configuration.__wbg_ptr, x, y);
    return Pair.__wrap(ret);
}

/**
* @param {FieldConfiguration} field_configuration
* @param {number} x
* @param {number} y
* @param {Float64Array} buffer
*/
export function compute_field_magnetostatic_direct_to_buffer(field_configuration, x, y, buffer) {
    _assertClass(field_configuration, FieldConfiguration);
    var ptr0 = passArrayF64ToWasm0(buffer, wasm.__wbindgen_malloc);
    var len0 = WASM_VECTOR_LEN;
    wasm.compute_field_magnetostatic_direct_to_buffer(field_configuration.__wbg_ptr, x, y, ptr0, len0, addHeapObject(buffer));
}

/**
* @param {FieldConfiguration} field_configuration
* @param {number} x
* @param {number} y
* @returns {Pair}
*/
export function compute_electric_field_dynamic(field_configuration, x, y) {
    _assertClass(field_configuration, FieldConfiguration);
    const ret = wasm.compute_electric_field_dynamic(field_configuration.__wbg_ptr, x, y);
    return Pair.__wrap(ret);
}

/**
* @param {FieldConfiguration} field_configuration
* @param {number} x
* @param {number} y
* @param {Float64Array} buffer
*/
export function compute_electric_field_dynamic_to_buffer(field_configuration, x, y, buffer) {
    _assertClass(field_configuration, FieldConfiguration);
    var ptr0 = passArrayF64ToWasm0(buffer, wasm.__wbindgen_malloc);
    var len0 = WASM_VECTOR_LEN;
    wasm.compute_electric_field_dynamic_to_buffer(field_configuration.__wbg_ptr, x, y, ptr0, len0, addHeapObject(buffer));
}

/**
* Generate a contour at a specified level of the electrostatic potential field
* @param {FieldConfiguration} field_configuration
* @param {Float64Array} levels
* @returns {any}
*/
export function generate_potential_contours_at_levels(field_configuration, levels) {
    _assertClass(field_configuration, FieldConfiguration);
    const ptr0 = passArrayF64ToWasm0(levels, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.generate_potential_contours_at_levels(field_configuration.__wbg_ptr, ptr0, len0);
    return takeObject(ret);
}

/**
* @param {FieldConfiguration} field_configuration
* @param {Float64Array} levels
* @returns {any}
*/
export function generate_potential_contours_and_arrow_positions_at_levels(field_configuration, levels) {
    _assertClass(field_configuration, FieldConfiguration);
    const ptr0 = passArrayF64ToWasm0(levels, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.generate_potential_contours_and_arrow_positions_at_levels(field_configuration.__wbg_ptr, ptr0, len0);
    return takeObject(ret);
}

function handleError(f, args) {
    try {
        return f.apply(this, args);
    } catch (e) {
        wasm.__wbindgen_exn_store(addHeapObject(e));
    }
}

const ChargeFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_charge_free(ptr >>> 0));
/**
*/
export class Charge {

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        ChargeFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_charge_free(ptr);
    }
    /**
    * @param {number} x
    * @param {number} y
    * @param {number} charge
    */
    constructor(x, y, charge) {
        const ret = wasm.charge_new(x, y, charge);
        this.__wbg_ptr = ret >>> 0;
        return this;
    }
}

const FieldConfigurationFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_fieldconfiguration_free(ptr >>> 0));
/**
*/
export class FieldConfiguration {

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        FieldConfigurationFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_fieldconfiguration_free(ptr);
    }
    /**
    * @param {number} x_max
    * @param {number} y_max
    * @param {number} nx
    * @param {number} ny
    */
    constructor(x_max, y_max, nx, ny) {
        const ret = wasm.fieldconfiguration_new(x_max, y_max, nx, ny);
        this.__wbg_ptr = ret >>> 0;
        return this;
    }
    /**
    * @param {any} charges
    */
    set_charges(charges) {
        wasm.fieldconfiguration_set_charges(this.__wbg_ptr, addHeapObject(charges));
    }
    /**
    * @param {number} x
    * @param {number} y
    */
    set_uniform_field(x, y) {
        wasm.fieldconfiguration_set_uniform_field(this.__wbg_ptr, x, y);
    }
    /**
    */
    reset_fields() {
        wasm.fieldconfiguration_reset_fields(this.__wbg_ptr);
    }
    /**
    */
    make_cic_grid() {
        wasm.fieldconfiguration_make_cic_grid(this.__wbg_ptr);
    }
    /**
    */
    initialize_on_constraints() {
        wasm.fieldconfiguration_initialize_on_constraints(this.__wbg_ptr);
    }
    /**
    */
    ensure_initialized() {
        wasm.fieldconfiguration_ensure_initialized(this.__wbg_ptr);
    }
    /**
    * @param {number} x
    * @param {number} y
    * @returns {number}
    */
    evaluate_cic_grid_interpolated(x, y) {
        const ret = wasm.fieldconfiguration_evaluate_cic_grid_interpolated(this.__wbg_ptr, x, y);
        return ret;
    }
    /**
    * @param {number} delta_t
    */
    make_currents(delta_t) {
        wasm.fieldconfiguration_make_currents(this.__wbg_ptr, delta_t);
    }
    /**
    * Evolve the fields by one timestep
    * @param {number} delta_t
    */
    tick(delta_t) {
        wasm.fieldconfiguration_tick(this.__wbg_ptr, delta_t);
    }
}

const PairFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_pair_free(ptr >>> 0));
/**
*/
export class Pair {

    static __wrap(ptr) {
        ptr = ptr >>> 0;
        const obj = Object.create(Pair.prototype);
        obj.__wbg_ptr = ptr;
        PairFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }

    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        PairFinalization.unregister(this);
        return ptr;
    }

    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_pair_free(ptr);
    }
    /**
    * @returns {number}
    */
    get u() {
        const ret = wasm.__wbg_get_pair_u(this.__wbg_ptr);
        return ret;
    }
    /**
    * @param {number} arg0
    */
    set u(arg0) {
        wasm.__wbg_set_pair_u(this.__wbg_ptr, arg0);
    }
    /**
    * @returns {number}
    */
    get v() {
        const ret = wasm.__wbg_get_pair_v(this.__wbg_ptr);
        return ret;
    }
    /**
    * @param {number} arg0
    */
    set v(arg0) {
        wasm.__wbg_set_pair_v(this.__wbg_ptr, arg0);
    }
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);

            } catch (e) {
                if (module.headers.get('Content-Type') != 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else {
                    throw e;
                }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);

    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };

        } else {
            return instance;
        }
    }
}

function __wbg_get_imports() {
    const imports = {};
    imports.wbg = {};
    imports.wbg.__wbindgen_string_new = function(arg0, arg1) {
        const ret = getStringFromWasm0(arg0, arg1);
        return addHeapObject(ret);
    };
    imports.wbg.__wbindgen_object_drop_ref = function(arg0) {
        takeObject(arg0);
    };
    imports.wbg.__wbindgen_copy_to_typed_array = function(arg0, arg1, arg2) {
        new Uint8Array(getObject(arg2).buffer, getObject(arg2).byteOffset, getObject(arg2).byteLength).set(getArrayU8FromWasm0(arg0, arg1));
    };
    imports.wbg.__wbindgen_error_new = function(arg0, arg1) {
        const ret = new Error(getStringFromWasm0(arg0, arg1));
        return addHeapObject(ret);
    };
    imports.wbg.__wbindgen_is_object = function(arg0) {
        const val = getObject(arg0);
        const ret = typeof(val) === 'object' && val !== null;
        return ret;
    };
    imports.wbg.__wbindgen_is_undefined = function(arg0) {
        const ret = getObject(arg0) === undefined;
        return ret;
    };
    imports.wbg.__wbindgen_in = function(arg0, arg1) {
        const ret = getObject(arg0) in getObject(arg1);
        return ret;
    };
    imports.wbg.__wbindgen_number_get = function(arg0, arg1) {
        const obj = getObject(arg1);
        const ret = typeof(obj) === 'number' ? obj : undefined;
        getFloat64Memory0()[arg0 / 8 + 1] = isLikeNone(ret) ? 0 : ret;
        getInt32Memory0()[arg0 / 4 + 0] = !isLikeNone(ret);
    };
    imports.wbg.__wbg_log_5bb5f88f245d7762 = function(arg0) {
        console.log(getObject(arg0));
    };
    imports.wbg.__wbg_log_1746d5c75ec89963 = function(arg0, arg1) {
        console.log(getObject(arg0), getObject(arg1));
    };
    imports.wbg.__wbg_new_abda76e883ba8a5f = function() {
        const ret = new Error();
        return addHeapObject(ret);
    };
    imports.wbg.__wbg_stack_658279fe44541cf6 = function(arg0, arg1) {
        const ret = getObject(arg1).stack;
        const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        getInt32Memory0()[arg0 / 4 + 1] = len1;
        getInt32Memory0()[arg0 / 4 + 0] = ptr1;
    };
    imports.wbg.__wbg_error_f851667af71bcfc6 = function(arg0, arg1) {
        let deferred0_0;
        let deferred0_1;
        try {
            deferred0_0 = arg0;
            deferred0_1 = arg1;
            console.error(getStringFromWasm0(arg0, arg1));
        } finally {
            wasm.__wbindgen_free(deferred0_0, deferred0_1, 1);
        }
    };
    imports.wbg.__wbindgen_jsval_loose_eq = function(arg0, arg1) {
        const ret = getObject(arg0) == getObject(arg1);
        return ret;
    };
    imports.wbg.__wbindgen_boolean_get = function(arg0) {
        const v = getObject(arg0);
        const ret = typeof(v) === 'boolean' ? (v ? 1 : 0) : 2;
        return ret;
    };
    imports.wbg.__wbindgen_string_get = function(arg0, arg1) {
        const obj = getObject(arg1);
        const ret = typeof(obj) === 'string' ? obj : undefined;
        var ptr1 = isLikeNone(ret) ? 0 : passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len1 = WASM_VECTOR_LEN;
        getInt32Memory0()[arg0 / 4 + 1] = len1;
        getInt32Memory0()[arg0 / 4 + 0] = ptr1;
    };
    imports.wbg.__wbindgen_number_new = function(arg0) {
        const ret = arg0;
        return addHeapObject(ret);
    };
    imports.wbg.__wbindgen_object_clone_ref = function(arg0) {
        const ret = getObject(arg0);
        return addHeapObject(ret);
    };
    imports.wbg.__wbg_getwithrefkey_15c62c2b8546208d = function(arg0, arg1) {
        const ret = getObject(arg0)[getObject(arg1)];
        return addHeapObject(ret);
    };
    imports.wbg.__wbg_get_bd8e338fbd5f5cc8 = function(arg0, arg1) {
        const ret = getObject(arg0)[arg1 >>> 0];
        return addHeapObject(ret);
    };
    imports.wbg.__wbg_length_cd7af8117672b8b8 = function(arg0) {
        const ret = getObject(arg0).length;
        return ret;
    };
    imports.wbg.__wbg_new_16b304a2cfa7ff4a = function() {
        const ret = new Array();
        return addHeapObject(ret);
    };
    imports.wbg.__wbindgen_is_function = function(arg0) {
        const ret = typeof(getObject(arg0)) === 'function';
        return ret;
    };
    imports.wbg.__wbg_next_40fc327bfc8770e6 = function(arg0) {
        const ret = getObject(arg0).next;
        return addHeapObject(ret);
    };
    imports.wbg.__wbg_next_196c84450b364254 = function() { return handleError(function (arg0) {
        const ret = getObject(arg0).next();
        return addHeapObject(ret);
    }, arguments) };
    imports.wbg.__wbg_done_298b57d23c0fc80c = function(arg0) {
        const ret = getObject(arg0).done;
        return ret;
    };
    imports.wbg.__wbg_value_d93c65011f51a456 = function(arg0) {
        const ret = getObject(arg0).value;
        return addHeapObject(ret);
    };
    imports.wbg.__wbg_iterator_2cee6dadfd956dfa = function() {
        const ret = Symbol.iterator;
        return addHeapObject(ret);
    };
    imports.wbg.__wbg_get_e3c254076557e348 = function() { return handleError(function (arg0, arg1) {
        const ret = Reflect.get(getObject(arg0), getObject(arg1));
        return addHeapObject(ret);
    }, arguments) };
    imports.wbg.__wbg_call_27c0f87801dedf93 = function() { return handleError(function (arg0, arg1) {
        const ret = getObject(arg0).call(getObject(arg1));
        return addHeapObject(ret);
    }, arguments) };
    imports.wbg.__wbg_set_d4638f722068f043 = function(arg0, arg1, arg2) {
        getObject(arg0)[arg1 >>> 0] = takeObject(arg2);
    };
    imports.wbg.__wbg_isArray_2ab64d95e09ea0ae = function(arg0) {
        const ret = Array.isArray(getObject(arg0));
        return ret;
    };
    imports.wbg.__wbg_instanceof_ArrayBuffer_836825be07d4c9d2 = function(arg0) {
        let result;
        try {
            result = getObject(arg0) instanceof ArrayBuffer;
        } catch (_) {
            result = false;
        }
        const ret = result;
        return ret;
    };
    imports.wbg.__wbg_buffer_12d079cc21e14bdb = function(arg0) {
        const ret = getObject(arg0).buffer;
        return addHeapObject(ret);
    };
    imports.wbg.__wbg_new_63b92bc8671ed464 = function(arg0) {
        const ret = new Uint8Array(getObject(arg0));
        return addHeapObject(ret);
    };
    imports.wbg.__wbg_set_a47bac70306a19a7 = function(arg0, arg1, arg2) {
        getObject(arg0).set(getObject(arg1), arg2 >>> 0);
    };
    imports.wbg.__wbg_length_c20a40f15020d68a = function(arg0) {
        const ret = getObject(arg0).length;
        return ret;
    };
    imports.wbg.__wbg_instanceof_Uint8Array_2b3bbecd033d19f6 = function(arg0) {
        let result;
        try {
            result = getObject(arg0) instanceof Uint8Array;
        } catch (_) {
            result = false;
        }
        const ret = result;
        return ret;
    };
    imports.wbg.__wbindgen_debug_string = function(arg0, arg1) {
        const ret = debugString(getObject(arg1));
        const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        getInt32Memory0()[arg0 / 4 + 1] = len1;
        getInt32Memory0()[arg0 / 4 + 0] = ptr1;
    };
    imports.wbg.__wbindgen_throw = function(arg0, arg1) {
        throw new Error(getStringFromWasm0(arg0, arg1));
    };
    imports.wbg.__wbindgen_memory = function() {
        const ret = wasm.memory;
        return addHeapObject(ret);
    };

    return imports;
}

function __wbg_init_memory(imports, maybe_memory) {

}

function __wbg_finalize_init(instance, module) {
    wasm = instance.exports;
    __wbg_init.__wbindgen_wasm_module = module;
    cachedFloat64Memory0 = null;
    cachedInt32Memory0 = null;
    cachedUint8Memory0 = null;


    return wasm;
}

function initSync(module) {
    if (wasm !== undefined) return wasm;

    const imports = __wbg_get_imports();

    __wbg_init_memory(imports);

    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }

    const instance = new WebAssembly.Instance(module, imports);

    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(input) {
    if (wasm !== undefined) return wasm;

    if (typeof input === 'undefined') {
        input = new URL('maxwell_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof input === 'string' || (typeof Request === 'function' && input instanceof Request) || (typeof URL === 'function' && input instanceof URL)) {
        input = fetch(input);
    }

    __wbg_init_memory(imports);

    const { instance, module } = await __wbg_load(await input, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync }
export default __wbg_init;
