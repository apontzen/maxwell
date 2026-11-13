// Playwright configuration for the MEME regression tests. See tests/README.md.
import { defineConfig, devices } from '@playwright/test';

const PORT = Number(process.env.MEME_TEST_PORT || 8123);

export default defineConfig({
    testDir: './tests',
    globalSetup: './tests/global-setup.js',
    fullyParallel: true,
    forbidOnly: !!process.env.CI,
    retries: 0,
    reporter: process.env.CI ? [['list'], ['html', { open: 'never' }]] : [['list']],

    // Reference images live next to the tests. They are platform-specific because
    // canvas anti-aliasing can differ between operating systems.
    snapshotPathTemplate: '{testDir}/__snapshots__/{testFileName}/{arg}-{platform}{ext}',

    expect: {
        toMatchSnapshot: {
            // `threshold` is the per-pixel colour tolerance (0-1) and `maxDiffPixels`
            // the number of pixels allowed to exceed it. Rendering is deterministic in
            // a fixed environment so these are tight: a single missing arrowhead or
            // a fieldline shifted by a pixel will fail.
            threshold: 0.1,
            maxDiffPixels: 10,
        },
    },

    use: {
        baseURL: `http://127.0.0.1:${PORT}`,
        trace: 'off',
        // Page screenshots stall while the fake clock is paused; canvas diffs are attached instead
        screenshot: 'off',
    },

    projects: [
        {
            name: 'chromium',
            use: {
                ...devices['Desktop Chrome'],
                viewport: { width: 900, height: 700 },
                deviceScaleFactor: 1,
            },
        },
    ],

    webServer: {
        command: `node tests/server.mjs ${PORT}`,
        url: `http://127.0.0.1:${PORT}/tests/harness.html`,
        reuseExistingServer: !process.env.CI,
        stdout: 'ignore',
    },
});
