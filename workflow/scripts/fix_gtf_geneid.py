#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

from utils import parse_gtf_attr, attd_to_str, gtfrow_to_str
from utils import GTFRow, GTFError
from utils import iterfile


def fix_gtf_geneid(args):


comments = []
rows = []

with iterfile(filepath) as fh:
    lines = (l.strip('\n') for l in fh)
    for i,l in enumerate(lines):
        if not l:
            continue
        if l.startswith('#'):
            comments.append((i,l))
        else:
            f = l.split('\t')
            f[-1] = parse_gtf_attr(f[-1])
            rows.append((i, GTFRow(*f)))

genes = defaultdict(list)
for lnum, r in rows:
    if 'gene_id' in r.attd:
        genes[r.attd['gene_id']].append((lnum, r))




with open('test_EBV.txt', 'w') as outh:
    for g, flist in genes.items():
        print(f'\nGENE: {g}', file=outh)
        for lnum, r in genes[g]:
            print(f'    {gtfrow_to_str(r)}', file=outh)












    """ Check that transcript_id exists for all """
for lnum,r in rows:
        if 'transcript_id' not in r.attd:
            raise GTFError('missing transcript_id:\n%s' % gtfrow_to_str(r))
    
    """ Mapping transcript_id to gene_id """
tid_to_gid = {}
for lnum, r in rows:
        if 'gene_id' in r.attd:
            if r.attd['transcript_id'] in tid_to_gid:
                if tid_to_gid[r.attd['transcript_id']] != r.attd['gene_id']:
                    raise GTFError(f'mismatch in line {lnum}:\n    {tid_to_gid[r.attd["transcript_id"]]} was seen, found {r.attd["gene_id"]}')
            tid_to_gid[r.attd['transcript_id']] = r.attd['gene_id']
    
    print('%d rows missing gene_id' % sum('gene_id' not in r.attd for r in rows), file=sys.stderr)
    
    """ Fix missing gene_id based on transcript_id"""
    for r in rows:
        if 'gene_id' not in r.attd:
            if r.attd['transcript_id'] in tid_to_gid:
                r.attd['gene_id'] = tid_to_gid[r.attd['transcript_id']]
            else:
                print('transcript_id not found: %s' % r.attd['transcript_id'])
    
    print('%d rows missing gene_id' % sum('gene_id' not in r.attd for r in rows), file=sys.stderr)
    
    for r in rows:
        r.attd.move_to_end('gene_id', last=False)
        r.attd.move_to_end('transcript_id', last=False)
        print('\t'.join(r['cols'] + [attd_to_str(r.attd)]), file=args.outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert GenBank/RefSeq GTF')
    parser.add_argument('--newref', type=str)
    parser.add_argument('infile')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args()                        
    fix_gtf_geneid(args)
