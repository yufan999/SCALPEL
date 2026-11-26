#!/usr/bin/env python3
"""Exercise the SCALPEL tagger and markdup on synthetic reads."""
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

TAGGER = Path(__file__).resolve().parent.parent / 'src/combined_cell_umi_tag.py'
SAMTOOLS = shutil.which('samtools')
BAM2BED = shutil.which('bam2bed')


def run(args, **kwargs):
    return subprocess.run([str(x) for x in args], check=True, text=True, capture_output=True, **kwargs)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        if not TAGGER.is_file():
            print(f'SKIP tagger/markdup: {TAGGER} not present')
            return
        if not SAMTOOLS:
            print('SKIP synthetic BAM: samtools unavailable')
            return
        version_output = subprocess.run([SAMTOOLS, '--version'], text=True, capture_output=True).stdout
        version_line = version_output.splitlines()[0] if version_output.splitlines() else ''
        version = re.match(r'^samtools\s+(\d+)\.(\d+)(?:\D|$)', version_line)
        if not version:
            print('SKIP synthetic BAM: unable to parse samtools version')
            return
        major, minor = map(int, version.groups())
        if (major, minor) < (1, 20):
            print(f'SKIP synthetic BAM: samtools {major}.{minor} found, 1.20 or newer required')
            return
        # Single-end flag 0 exercises markdup without paired-read fixmate/MC/ms requirements.
        header = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n'
        def rec(name, pos, tags):
            return f'{name}\t0\tchr1\t{pos}\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\t' + '\t'.join(tags) + '\n'
        a, b, u, v = 'CB:Z:AAAA-1', 'CB:Z:BBBB-1', 'UB:Z:TTTT', 'UB:Z:GGGG'
        sam = header + rec('a', 100, [a, u]) + rec('b', 100, [b, u]) + rec('dup', 100, [a, u]) + rec('other', 100, [a, v]) + rec('no_ub', 200, [a]) + rec('no_cb', 300, [u])
        source = d / 'input.sam'
        source.write_text(sam)
        bam = d / 'input.bam'
        run([SAMTOOLS, 'view', '-b', '-o', bam, source])
        run([SAMTOOLS, 'index', bam])
        whitelist = d / 'barcodes'
        whitelist.write_text('AAAA-1\nBBBB-1\n')
        # Binary pipeline mirrors whitelist filtering and coordinate sorting.
        subprocess.run(f'"{SAMTOOLS}" view -b -D CB:"{whitelist}" --keep-tag CB,UB "{bam}" | "{SAMTOOLS}" sort -o "{d}/tmp.bam"', shell=True, executable='/bin/bash', check=True)
        tagged = run([SAMTOOLS, 'view', '-h', d / 'tmp.bam'])
        transformed = run([os.sys.executable, TAGGER], input=tagged.stdout)
        assert json.loads(transformed.stderr) == {'CB_UB_tagged_reads': 4, 'reads_without_CB_or_UB': 1}
        # Missing CB is filtered by whitelist before the helper; test its count directly too.
        missing = run([os.sys.executable, TAGGER], input=header + rec('no_cb', 300, [u]) + rec('no_ub', 200, [a]))
        assert json.loads(missing.stderr)['reads_without_CB_or_UB'] == 2
        (d / 'tagged.sam').write_text(transformed.stdout)
        run([SAMTOOLS, 'view', '-b', '-o', d / 'tagged.bam', d / 'tagged.sam'])
        run([SAMTOOLS, 'markdup', d / 'tagged.bam', d / 'dedup.bam', '-r', '--barcode-tag', 'XU'])
        run([SAMTOOLS, 'view', '-b', '-x', 'XU', '-o', d / 'final.bam', d / 'dedup.bam'])
        records = [x.split('\t') for x in run([SAMTOOLS, 'view', d / 'final.bam']).stdout.splitlines()]
        assert len(records) == 3
        assert {(x[11], x[12]) for x in records} == {(a, u), (b, u), (a, v)}
        assert all([z[:2] for z in x[11:]] == ['CB', 'UB'] for x in records)
        run([SAMTOOLS, 'markdup', d / 'tmp.bam', d / 'upstream.bam', '-r', '--barcode-tag', 'CB', '--barcode-tag', 'UB'])
        upstream = run([SAMTOOLS, 'view', d / 'upstream.bam']).stdout.splitlines()
        assert sum('UB:Z:TTTT' in x for x in upstream) == 1
        if BAM2BED:
            with (d / 'final.bam').open('rb') as src:
                bed = subprocess.run([str(BAM2BED), '--all-reads', '--split', '--do-not-sort'], stdin=src, text=True, capture_output=True, check=True, env={**os.environ, 'PATH': str(Path(SAMTOOLS).parent) + os.pathsep + os.environ.get('PATH', '')}).stdout.splitlines()
            assert {(x.split('\t')[13], x.split('\t')[14]) for x in bed} == {(a, u), (b, u), (a, v)}
        else:
            print('SKIP bam2bed field test: bam2bed unavailable')
        print('PASS: tagger, synthetic dedup, tags and bam2bed fields')


if __name__ == '__main__':
    main()
