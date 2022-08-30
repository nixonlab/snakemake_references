#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

from utils import parse_gtf_attr, attd_to_str, GTFError


def fix_gtf_geneid(args):
    lines = (l.strip().split('\t') for l in args.infile if not l.startswith('#'))
    rows = [{'cols': l[:-1], 'attd': parse_gtf_attr(l[-1])} for l in lines ]
    
    """ Check that transcript_id exists for all """
    for r in rows:
        if 'transcript_id' not in r['attd']:
            raise GTFError('missing transcript_id:\n%s' % '\t'.join(r['cols'] + [attd_to_str(r['attd'])]))
    
    """ Mapping transcript_id to gene_id """
    tid_to_gid = {}
    for r in rows:
        if 'gene_id' in r['attd']:
            if r['attd']['transcript_id'] in tid_to_gid:
                if tid_to_gid[r['attd']['transcript_id']] != r['attd']['gene_id']:
                    raise GTFError('mismatch: %s %s' % (tid_to_gid[r['attd']['transcript_id']], r['attd']['gene_id']))
            tid_to_gid[r['attd']['transcript_id']] = r['attd']['gene_id']
    
    print('%d rows missing gene_id' % sum('gene_id' not in r['attd'] for r in rows), file=sys.stderr)
    
    """ Fix missing gene_id based on transcript_id"""
    for r in rows:
        if 'gene_id' not in r['attd']:
            if r['attd']['transcript_id'] in tid_to_gid:
                r['attd']['gene_id'] = tid_to_gid[r['attd']['transcript_id']]
            else:
                print('transcript_id not found: %s' % r['attd']['transcript_id'])
    
    print('%d rows missing gene_id' % sum('gene_id' not in r['attd'] for r in rows), file=sys.stderr)
    
    for r in rows:
        r['attd'].move_to_end('gene_id', last=False)
        r['attd'].move_to_end('transcript_id', last=False)
        print('\t'.join(r['cols'] + [attd_to_str(r['attd'])]), file=args.outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix missing "gene_id" attributes in GTF file.')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args()                        
    fix_gtf_geneid(args)
