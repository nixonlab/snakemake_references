#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import pandas as pd

from utils import GTFCOLS, convert_attstr

# MEMPROF = False
# if MEMPROF:
#     from utils import memory_usage

# import argparse
# args = argparse.Namespace(
#     infile='databases/annotations/gencode.v38/gencode.v38.ALL.annotation.gtf.gz',
#     out_prefix='databases/annotations/gencode.v38/testrun',
#     intfeats=['exon', 'gene']
# )

def get_metadata(args):
    intfeats = ['gene', 'transcript']
    if args.intfeats:
        intfeats = list(pd.Series(intfeats + args.intfeats)
                        .drop_duplicates()
                        )

    print(f'Feature types to process: {", ".join(intfeats)}.', file=sys.stderr)

    allfeats = pd.read_csv(
        args.infile,
        sep='\t',
        comment='#',
        names=GTFCOLS,
        converters={'start':int, 'end':int, 'attributes': convert_attstr},
    )
    print("Loaded GTF", file=sys.stderr)
    # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
    # if MEMPROF: print(memory_usage(), file=sys.stderr)

    # Subset gene, transcript, and other features
    feattypes = {}
    for ftype in intfeats:
        _feats = allfeats[allfeats['feature']==ftype]
        _feats = pd.concat(
            [
                _feats.drop(columns='attributes'),
                _feats['attributes'].apply(pd.Series)
            ],
            axis=1
        )
        feattypes[ftype] = _feats
        # print(f"Completed parsing {ftype} features", file=sys.stderr)
        # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
        # if MEMPROF: print(memory_usage(), file=sys.stderr)

    del allfeats
    # print("Deleted allfeats", file=sys.stderr)
    # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
    # if MEMPROF: print(memory_usage(), file=sys.stderr)

    # make output table for other feature types
    ofeats = [f for f in intfeats if f not in ['gene', 'transcript']]
    for ftype in ofeats:
        _outfn = f'{args.out_prefix}.metadata.{ftype}_features.txt.gz'
        _feats = feattypes.pop(ftype)
        if ftype == 'exon':
            _feats.set_index('exon_id').fillna('').to_csv(_outfn, sep='\t')
        else:
            _feats.fillna('').to_csv(_outfn, sep='\t', index=False)
        del _feats
        print(f"Wrote metadata for '{ftype}' to {_outfn}", file=sys.stderr)
        # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
        # if MEMPROF: print(memory_usage(), file=sys.stderr)

    # Gene features
    genefeats = feattypes.pop('gene')

    # metadata.gene_features
    _outfn = f'{args.out_prefix}.metadata.gene_features.txt.gz'
    genefeats.set_index('gene_id').fillna('').to_csv(_outfn, sep='\t')
    print(f"Wrote metadata for 'gene' to {_outfn}", file=sys.stderr)

    # gene ID -> gene name
    _outfn = f'{args.out_prefix}.metadata.gid_gname.txt.gz'
    gid_gname = genefeats[['gene_id', 'gene_name']].drop_duplicates()
    if gid_gname['gene_id'].duplicated().sum():
        sys.stderr.write('WARNING: gene_id is repeated in gid_gname\n')
    gid_gname.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'gid_gname' to {_outfn}", file=sys.stderr)
    del gid_gname

    # print(f"Deleted gid_gname", file=sys.stderr)
    # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
    # if MEMPROF: print(memory_usage(), file=sys.stderr)

    # gene ID -> gene type
    _outfn = f'{args.out_prefix}.metadata.gid_gtype.txt.gz'
    gid_gtype = genefeats[['gene_id', 'gene_type']].drop_duplicates()
    if gid_gtype['gene_id'].duplicated().sum():
        sys.stderr.write('WARNING: gene_id is repeated in gid_gtype\n')
    if gid_gtype['gene_type'].isnull().sum():
        sys.stderr.write('WARNING: gene_type is NULL in gid_gtype\n')
    gid_gtype.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'gid_gtype' to {_outfn}", file=sys.stderr)
    del gid_gtype

    # print(f"Deleted gid_gtype", file=sys.stderr)
    # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
    # if MEMPROF: print(memory_usage(), file=sys.stderr)

    # gene ID -> hgnc
    _outfn = f'{args.out_prefix}.metadata.gid_hgnc.txt.gz'
    gid_hgnc = genefeats[['gene_id', 'hgnc_id']].drop_duplicates()
    if gid_hgnc['gene_id'].duplicated().sum():
        sys.stderr.write('WARNING: gene_id is repeated in gid_hgnc\n')
    # if gid_hgnc['hgnc_id'].isnull().sum():
    #     sys.stderr.write('WARNING: hgnc_id is NULL in gid_hgnc\n')
    gid_hgnc.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'gid_hgnc' to {_outfn}", file=sys.stderr)
    del gid_hgnc

    # print(f"Deleted gid_hgnc", file=sys.stderr)
    # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
    # if MEMPROF: print(memory_usage(), file=sys.stderr)

    # Done with genefeats
    del genefeats

    # print(f"Deleted genefeats", file=sys.stderr)
    # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
    # if MEMPROF: print(memory_usage(), file=sys.stderr)

    # Transcript features
    _outfn = f'{args.out_prefix}.metadata.tx_features.txt.gz'
    txfeats = feattypes.pop('transcript')
    txfeats.set_index('transcript_id').fillna('').to_csv(_outfn, sep='\t')
    print(f"Wrote metadata for 'transcript' to {_outfn}", file=sys.stderr)

    # transcript ID -> transcript name
    _outfn = f'{args.out_prefix}.metadata.tid_tname.txt.gz'
    tid_tname = txfeats[['transcript_id', 'transcript_name']].drop_duplicates()
    if tid_tname['transcript_id'].duplicated().sum():
        sys.stderr.write('WARNING: transcript_id is repeated in tid_tname\n')
    tid_tname.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'tid_tname' to {_outfn}", file=sys.stderr)
    del tid_tname

    # print(f"Deleted tid_tname", file=sys.stderr)
    # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
    # if MEMPROF: print(memory_usage(), file=sys.stderr)

    # transcript ID -> transcript type
    _outfn = f'{args.out_prefix}.metadata.tid_ttype.txt.gz'
    tid_ttype = txfeats[['transcript_id', 'transcript_type']].drop_duplicates()
    if tid_ttype['transcript_id'].duplicated().sum():
        sys.stderr.write('WARNING: transcript_id is repeated in tid_ttype\n')
    if tid_ttype['transcript_type'].isnull().sum():
        sys.stderr.write('WARNING: transcript_type is NULL in tid_ttype\n')
    tid_ttype.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'tid_ttype' to {_outfn}", file=sys.stderr)
    del tid_ttype

    # print(f"Deleted tid_ttype", file=sys.stderr)
    # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
    # if MEMPROF: print(memory_usage(), file=sys.stderr)

    # transcript ID -> hgnc
    _outfn = f'{args.out_prefix}.metadata.tid_hgnc.txt.gz'
    tid_hgnc = txfeats[['transcript_id', 'hgnc_id']].drop_duplicates()
    if tid_hgnc['transcript_id'].duplicated().sum():
        sys.stderr.write('WARNING: transcript_id is repeated in tid_hgnc\n')
    # if tid_hgnc['hgnc_id'].isnull().sum():
    #     sys.stderr.write('WARNING: hgnc_id is NULL in tid_hgnc\n')
    tid_hgnc.dropna().to_csv(_outfn, sep='\t', header=False, index=False)
    print(f"  Wrote 'tid_hgnc' to {_outfn}", file=sys.stderr)
    del tid_hgnc

    # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
    # if MEMPROF: print(memory_usage(), file=sys.stderr)

    # transcript_id to gene_id
    _outfn = f'{args.out_prefix}.metadata.tid_gid.txt.gz'
    tid_gid = txfeats[['transcript_id', 'gene_id']].drop_duplicates()
    if tid_gid['transcript_id'].duplicated().sum():
        sys.stderr.write('WARNING: transcript_id is repeated in tid_gid\n')
    if tid_gid['gene_id'].isnull().sum():
        sys.stderr.write('WARNING: gene_id is NULL in tid_gid\n')
    tid_gid.set_index('transcript_id').fillna('').to_csv(_outfn, sep='\t', header=False)
    print(f"  Wrote 'tid_gid' to {_outfn}", file=sys.stderr)

    # gene_id to transcript_id list
    _outfn = f'{args.out_prefix}.metadata.gid_tid.txt.gz'
    gid_tid = tid_gid.groupby('gene_id')['transcript_id'].apply(list)
    gid_tid.apply(lambda x:','.join(map(str,x))).to_csv(_outfn, sep='\t', header=False)
    print(f"  Wrote 'gid_tid' to {_outfn}", file=sys.stderr)

    # if MEMPROF: print("Calculating memory usage", file=sys.stderr)
    # if MEMPROF: print(memory_usage(), file=sys.stderr)


"""
def old_get_metadata(args):
    if args.gff:
        iter = (GFFRow(row) for row in args.infile if not row.startswith('#'))
    else:
        iter = (GTFRow(row) for row in args.infile if not row.startswith('#'))

    tid_gid = {}
    gid_tids = defaultdict(set)
    gid_gtype = {}
    gid_gname = {}
    tid_ttype = {}
    tid_tname = {}

    for i, g in enumerate(iter):
        if g.feature == 'gene':
            assert 'gene_id' in g.attd, f'gene_id not in "gene" feature, {i}\n{g}'
            _gid = g.attd['gene_id']
            if 'gene_name' in g.attd:
                if not check_setdefault(gid_gname, _gid, g.attd['gene_name']):
                    print(f'gene_name mismatch ({_gid})\nrow{i}: {g}',
                          file=sys.stderr
                    )
            if 'gene_type' in g.attd:
                if not check_setdefault(gid_gname, _gid, g.attd['gene_type']):
                    print(f'gene_type mismatch ({_gid})\nrow{i}: {g}',
                          file=sys.stderr
                    )
        elif g.feature == 'transcript':
            assert 'gene_id' in g.attd, f'gene_id not in "transcript" feature, {i}\n{g}'
            assert 'transcript_id' in g.attd, f'transcript_id not in "transcript" feature, {i}\n{g}'
            _gid = g.attd['gene_id']
            _tid = g.attd['transcript_id']
            check_setdefault(tid_gid, _tid, _gid, i, g)
            gid_tids[_gid].add(_tid)
            if 'transcript_name' in g.attd:
                if not check_setdefault(gid_gname, _gid, g.attd['transcript_name']):
                    print(f'transcript_name mismatch ({_gid})\nrow{i}: {g}',
                          file=sys.stderr
                    )
            if 'transcript_type' in g.attd:
                if not check_setdefault(gid_gname, _gid, g.attd['transcript_type']):
                    print(f'transcript_type mismatch ({_gid})\nrow{i}: {g}',
                          file=sys.stderr
                    )
            #
            # if 'transcript_type' in g.attd:
            #     check_setdefault(tid_ttype, g.attd['transcript_id'], g.attd['transcript_type'], i, g)
            # if 'transcript_name' in g.attd:
            #     check_setdefault(tid_tname, g.attd['transcript_id'], g.attd['transcript_name'], i, g)
        #if 'gene_id' in g.attd:
        #    if 'gene_type' in g.attd:
        #        check_setdefault(gid_gtype, g.attd['gene_id'], g.attd['gene_type'], i, g)
        #    if 'gene_name' in g.attd:
        #        check_setdefault(gid_gname, g.attd['gene_id'], g.attd['gene_name'], i, g)
        #if 'transcript_id' in g.attd:
        #    if 'transcript_type' in g.attd:
        #         check_setdefault(tid_ttype, g.attd['transcript_id'], g.attd['transcript_type'], i, g)
        #     if 'transcript_name' in g.attd:
        #         check_setdefault(tid_tname, g.attd['transcript_id'], g.attd['transcript_name'], i, g)

    with gzip.open(f'{args.out_prefix}.metadata.tid_gid.txt.gz', 'wb') as outh:
        for tid,gid in tid_gid.items():
            outh.write(f'{tid}\t{gid}\n'.encode())

    with gzip.open(f'{args.out_prefix}.metadata.tid_ttype.txt.gz', 'wb') as outh:
        for tid,ttype in tid_ttype.items():
            outh.write(f'{tid}\t{ttype}\n'.encode())

    with gzip.open(f'{args.out_prefix}.metadata.tid_tname.txt.gz', 'wb') as outh:
        for tid,tname in tid_tname.items():
            outh.write(f'{tid}\t{tname}\n'.encode())

    with gzip.open(f'{args.out_prefix}.metadata.gid_gtype.txt.gz', 'wb') as outh:
        for gid,gtype in gid_gtype.items():
            outh.write(f'{gid}\t{gtype}\n'.encode())

    with gzip.open(f'{args.out_prefix}.metadata.gid_gname.txt.gz', 'wb') as outh:
        for gid,gname in gid_gname.items():
            outh.write(f'{gid}\t{gname}\n'.encode())

    with gzip.open(f'{args.out_prefix}.metadata.gid_tid.txt.gz', 'wb') as outh:
        for gid, tids in gid_tids.items():
            outh.write(f'{gid}\t{",".join(sorted(tids))}\n'.encode())

    if args.tid_HGNC is None:
        return

    gid_hgnc = defaultdict(Counter)
    hgnc_gid = defaultdict(Counter)
    with gzip.open(args.tid_HGNC, 'rt') as fh:
        _iter = (_.strip('\n').split('\t') for _ in fh)
        for tid, hgnc_name, hgnc_id in _iter:
            try:
                gid = tid_gid[tid]
                gid_hgnc[gid][(hgnc_name, hgnc_id)] += 1
                hgnc_gid[hgnc_name][gid] += 1
            except KeyError as e:
                msg = f'WARNING: transcript ID {id} from tid_HGNC was not in GTF'
                print(msg, file = sys.stderr)

    with gzip.open(f'{args.out_prefix}.metadata.gid_HGNC.txt.gz', 'wb') as outh:
        for gid, vals in gid_hgnc.items():
            if len(vals) == 1:
                (hgnc_name, hgnc_id), n = vals.popitem()
                outh.write(f'{gid}\t{hgnc_name}\t{hgnc_id}\n'.encode())
            else:
                _s = '; '.join(f'{k[0]}/{k[1]} ({ct})' for k, ct in vals.most_common())
                msg = f'WARNING: gene ID {gid} matched to multiple names:\t{_s}'
                print(msg, file = sys.stderr)
                (hgnc_name, hgnc_id), n = vals.most_common()[0]
                outh.write(f'{gid}\t{hgnc_name}\t{hgnc_id}\n'.encode())

    with gzip.open(f'{args.out_prefix}.metadata.HGNC_gid.txt.gz', 'wb') as outh:
        for hgnc_name, vals in hgnc_gid.items():
            if len(vals) == 1:
                gid, n = vals.popitem()
                outh.write(f'{hgnc_name}\t{gid}\n'.encode())
            else:
                _s = '; '.join(f'{k} ({ct})' for k, ct in vals.most_common())
                msg = f'WARNING: {hgnc_name} matched multiple gene IDs:\t{_s}'
                print(msg, file = sys.stderr)
                # csv method
                # gid_csv = ','.join(k for k,ct in vals.most_common())
                # outh.write(f'{hgnc_name}\t{gid_csv}\n'.encode())
                # one per line
                for gid,ct in vals.items():
                    outh.write(f'{hgnc_name}\t{gid}\n'.encode())
"""


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get metadata from GTF')
    parser.add_argument('infile',
                        help="Input GTF."
                        )
    parser.add_argument('out_prefix',
                        help = "Output file prefix."
                        )
    parser.add_argument('--intfeats', action='append',
                        help = '''Create metadata tables for these feature 
                                  types (can specify multiple times). Tables 
                                  for 'gene' and 'transcript' 
                                  features are always enabled, other feature
                                  types may include: ['exon', 'CDS', 'UTR', 
                                  'start_codon', 'stop_codon', 
                                  'Selenocysteine'].
                                  ''')
    args = parser.parse_args()
    get_metadata(args)

args = parser.parse_args()