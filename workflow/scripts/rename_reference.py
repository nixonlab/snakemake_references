#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys

from utils import GTFRow, GFFRow

def rename_reference(args):
    GXFRow = GFFRow if args.gff else GTFRow
    for row in args.infile:
        if row.startswith('#'):
            print(row.strip('\n'), file=args.outfile)
            continue
        gtfrow = GXFRow(row)
        gtfrow.attd['orig_ref'] = gtfrow.chrom
        gtfrow.chrom = args.newref
        print(gtfrow, file=args.outfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Rename reference')
    parser.add_argument('--gff', action='store_true')
    parser.add_argument('newref', type=str)
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args()
    rename_reference(args)
