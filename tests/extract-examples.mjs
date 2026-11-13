#!/usr/bin/env node
// Regenerate tests/fixtures/examples.json from the <meme-embed> elements that
// appear in the lecture notes (apontzen/du_electromagnetism) and in embed_test.html.
//
// Usage:
//   node tests/extract-examples.mjs [path/to/du_electromagnetism]
//
// If no path is given, ../du_electromagnetism (a sibling checkout) is tried, then
// ../.. (i.e. this repo checked out as the MEME submodule of du_electromagnetism).
// Embeds inside HTML comments are skipped, since they are not shown to students.

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(here, '..');

/** Locate a checkout of du_electromagnetism, or return null if there isn't one. */
export function findLectureRepo(explicitPath = null) {
    const candidates = explicitPath
        ? [explicitPath]
        : [path.resolve(repoRoot, '../du_electromagnetism'), path.resolve(repoRoot, '..')];
    for (const c of candidates) {
        if (fs.existsSync(path.join(c, 'lectures')) && fs.existsSync(path.join(c, '_config.yml')))
            return path.resolve(c);
    }
    return null;
}

function decodeEntities(s) {
    return s.replace(/&quot;/g, '"').replace(/&#39;/g, "'").replace(/&apos;/g, "'")
        .replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&amp;/g, '&');
}

function lineOf(text, index) {
    return text.slice(0, index).split('\n').length;
}

// Return [start, end) ranges of HTML comments so that embeds inside them can be skipped
function commentRanges(text) {
    const ranges = [];
    const re = /<!--[\s\S]*?-->/g;
    let m;
    while ((m = re.exec(text)) !== null) ranges.push([m.index, m.index + m[0].length]);
    return ranges;
}

function parseAttributes(tagBody) {
    const attrs = {};
    const re = /([a-zA-Z_:][-a-zA-Z0-9_:.]*)\s*(?:=\s*(?:"([^"]*)"|'([^']*)'|([^\s"'>]+)))?/g;
    let m;
    while ((m = re.exec(tagBody)) !== null) {
        const value = m[2] ?? m[3] ?? m[4] ?? '';
        attrs[m[1]] = decodeEntities(value);
    }
    return attrs;
}

function extract(file, relName, idPrefix) {
    const text = fs.readFileSync(file, 'utf8');
    const comments = commentRanges(text);
    const re = /<meme-embed\b([\s\S]*?)>/g;
    const results = [];
    let m, n = 0;
    while ((m = re.exec(text)) !== null) {
        if (comments.some(([a, b]) => m.index >= a && m.index < b)) continue;
        n++;
        const attributes = parseAttributes(m[1]);
        // Validate that the state really is JSON, and store a parsed copy for the tests' convenience
        const state = JSON.parse(attributes.meme);
        results.push({
            id: `${idPrefix}-${n}`,
            source: `${relName}:${lineOf(text, m.index)}`,
            attributes,
            state,
        });
    }
    return results;
}

/** Extract all examples from the lecture notes in `lectureRepo` and from embed_test.html */
export function extractExamples(lectureRepo) {
    const lectureDir = path.join(lectureRepo, 'lectures');
    const examples = [];
    for (const f of fs.readdirSync(lectureDir).filter(f => f.endsWith('.md')).sort()) {
        examples.push(...extract(path.join(lectureDir, f), `du_electromagnetism/lectures/${f}`, path.basename(f, '.md')));
    }
    examples.push(...extract(path.join(repoRoot, 'embed_test.html'), 'embed_test.html', 'embed_test'));
    return examples;
}

export const fixturePath = path.join(here, 'fixtures', 'examples.json');

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
    const lectureRepo = findLectureRepo(process.argv[2]);
    if (!lectureRepo) {
        console.error('Could not find a du_electromagnetism checkout; pass its path as an argument.');
        process.exit(1);
    }
    const examples = extractExamples(lectureRepo);
    fs.writeFileSync(fixturePath, JSON.stringify(examples, null, 2) + '\n');
    console.log(`Wrote ${examples.length} examples from ${lectureRepo} to ${path.relative(process.cwd(), fixturePath)}`);
}
