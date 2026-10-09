#!/usr/bin/env python3
"""Replace <!--SCATTER:key--> placeholders in a draft report with the
corresponding CNVkit scatter plot (from CREATE_CNVKIT_SCATTER's manifest-driven
output), producing the final report markdown."""

import argparse
import csv
import os
import re
import sys

PLACEHOLDER_RE = re.compile(r'<!--SCATTER:([^>]+)-->')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report',      required=True, help='Draft report markdown (from GENERATE_REPORT)')
    p.add_argument('--manifest',    required=True, help='Scatter manifest TSV (from GENERATE_REPORT)')
    p.add_argument('--scatter-dir', required=True, help='Directory of {key}.png scatter plots (from CREATE_CNVKIT_SCATTER)')
    p.add_argument('--out',         required=True, help='Final report markdown, with images embedded')
    return p.parse_args()


def load_manifest_keys(manifest_path):
    keys = []
    with open(manifest_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row.get('key'):
                keys.append(row['key'])
    return keys


def main():
    args = parse_args()

    manifest_keys = set(load_manifest_keys(args.manifest))

    with open(args.report) as f:
        report_text = f.read()

    placeholder_keys = set(PLACEHOLDER_RE.findall(report_text))
    missing_from_manifest = placeholder_keys - manifest_keys
    if missing_from_manifest:
        print(f"WARNING: {len(missing_from_manifest)} placeholder(s) in report have no manifest entry: "
              f"{sorted(missing_from_manifest)}", file=sys.stderr)

    n_embedded, n_unavailable = 0, 0

    def _replace(match):
        nonlocal n_embedded, n_unavailable
        key = match.group(1)
        png_path = os.path.join(args.scatter_dir, f"{key}.png")
        if os.path.isfile(png_path):
            n_embedded += 1
            return f"![{key} scatter](scatter/{key}.png)"
        n_unavailable += 1
        print(f"WARNING: scatter plot missing for key '{key}' (expected {png_path})", file=sys.stderr)
        return f"_(scatter plot unavailable for {key})_"

    final_text = PLACEHOLDER_RE.sub(_replace, report_text)

    with open(args.out, 'w') as f:
        f.write(final_text)

    print(f"Final report written to {args.out} ({n_embedded} scatter plots embedded, "
          f"{n_unavailable} unavailable)", file=sys.stderr)


if __name__ == '__main__':
    main()
