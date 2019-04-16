#!/usr/bin/env python3
"""
Copyright 2019 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Compare-annotations/

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not,
see <https://www.gnu.org/licenses/>.
"""

import argparse
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import itertools
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Compare GenBank annotations',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('annotation_1', type=str,
                        help='First annotated genome in Genbank format')
    parser.add_argument('annotation_2', type=str,
                        help='Second annotated genome in Genbank format')

    parser.add_argument('--match_identity_threshold', type=float, default=0.7,
                        help='Two genes must have at least this identity to be considerd the same (0.0 to 1.0)')
    parser.add_argument('--allowed_skipped_genes', type=int, default=10,
                        help='This many missing genes are allowed when aligning the annotations')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()

    # Load in the CDS features from the two assemblies.
    old = SeqIO.parse(args.annotation_1, 'genbank')
    new = SeqIO.parse(args.annotation_2, 'genbank')
    old_record = next(old)
    new_record = next(new)
    old_features, new_features = [], []
    for f in old_record.features:
        if f.type == 'CDS':
            old_features.append(f)
    for f in new_record.features:
        if f.type == 'CDS':
            new_features.append(f)

    print('Features in old assembly:', len(old_features))
    print('Features in new assembly:', len(new_features))

    # Align the features to each other.
    offsets = sorted(list(itertools.product(range(args.allowed_skipped_genes),
                                            range(args.allowed_skipped_genes))),
                     key=lambda x: x[0]+x[1])
    aligned_features = []
    old_i, new_i = 0, 0
    while True:
        if old_i >= len(old_features) and new_i >= len(new_features):
            break

        for old_offset, new_offset in offsets:
            try:
                old_feature = old_features[old_i + old_offset]
            except IndexError:
                old_feature = None
            try:
                new_feature = new_features[new_i + new_offset]
            except IndexError:
                new_feature = None

            match, identity, length_diff = compare_features(old_feature, new_feature, old_record, new_record, args.match_identity_threshold)
            if match:
                for j in range(old_offset):
                    print_in_old_not_new(old_features[old_i + j])
                for j in range(new_offset):
                    print_in_new_not_old(new_features[new_i + j])
                print_match(old_features[old_i + old_offset], new_features[new_i + new_offset], identity, length_diff)
                old_i += old_offset
                new_i += new_offset
                break
        else:
            sys.exit('\nERROR: Failed to find alignment')

        if old_feature is None and new_feature is None:
            break

        old_i += 1
        new_i += 1


def print_match(f1, f2, identity, length_diff):
    print('', flush=True)
    if identity == 1.0:
        print('Exact match')
    else:
        print('Inexact match (' + '%.2f' % (identity * 100.0) + '% ID, ', end='')
        if length_diff == 0:
            print('same length)')
        elif length_diff > 0:
            print('old seq longer)')
        elif length_diff < 0:
            print('new seq longer)')
    print('  old: ', end='')
    print_feature_one_line(f1)
    print('  new: ', end='')
    print_feature_one_line(f2)
    p1 = f1.qualifiers['product'][0].lower()
    p2 = f2.qualifiers['product'][0].lower()
    if 'hypothetical' in p1 and 'hypothetical' in p2:
        print('  still hypothetical')
    if 'hypothetical' in p1 and 'hypothetical' not in p2:
        print('  no longer hypothetical')
    if 'hypothetical' not in p1 and 'hypothetical' in p2:
        print('  became hypothetical')


def print_in_old_not_new(f):
    print('')
    print('In old but not in new:')
    print('  ', end='')
    print_feature_one_line(f)


def print_in_new_not_old(f):
    print('')
    print('In new but not in old:')
    print('  ', end='')
    print_feature_one_line(f)


def print_feature_one_line(f):
    f_str = f.qualifiers['product'][0]
    strand = '+' if f.location.strand == 1 else '-'
    f_str += ' (' + str(f.location.start) + '-' + str(f.location.end) + ' ' + strand + ', '
    f_str += str(f.location.end - f.location.start) + ' bp)'
    print(f_str)


def compare_features(f1, f2, r1, r2, match_identity_threshold):
    if f1 is None or f2 is None:
        return False
    s1 = f1.extract(r1).seq
    s2 = f2.extract(r2).seq
    score = pairwise2.align.globalms(s1, s2, 1, 0, 0, 0, score_only=True)
    identity = score / max(len(s1), len(s2))
    match = identity > match_identity_threshold
    length_diff = len(s1) - len(s2)
    return match, identity, length_diff


if __name__ == '__main__':
    main()
