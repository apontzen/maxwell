import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

export default function globalSetup() {
    const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
    const wasm = path.join(root, 'maxwell', 'out', 'maxwell_bg.wasm');
    if (!fs.existsSync(wasm)) {
        throw new Error(`${wasm} not found. Build the WebAssembly module first with \`npm run build:wasm\`.`);
    }
}
