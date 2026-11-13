// Checks that fixtures/examples.json is in sync with the lecture notes, when a
// du_electromagnetism checkout is available (either alongside this repository, or
// as its parent when this repository is the MEME submodule). Skipped otherwise.

import fs from 'node:fs';
import { test, expect } from '@playwright/test';
import { findLectureRepo, extractExamples, fixturePath } from './extract-examples.mjs';

test('fixtures/examples.json matches the <meme-embed>s in the lecture notes', () => {
    const lectureRepo = findLectureRepo(process.env.DU_ELECTROMAGNETISM || null);
    test.skip(!lectureRepo, 'no du_electromagnetism checkout found');
    const current = JSON.parse(fs.readFileSync(fixturePath, 'utf8'));
    expect(extractExamples(lectureRepo),
        'lecture notes have changed: run `npm run extract-examples` then `npm run test:update` ' +
        'and review the new/changed reference images').toEqual(current);
});
