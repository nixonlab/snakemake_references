#! /usr/bin/env python
# -*- coding: utf-8 -*-

import re
from collections import OrderedDict
import gzip

from typing import Optional


class GTFError(Exception):
    pass


class GTFRow(object):
    def __init__(self, row: Optional[str]):
        fields = row.strip('\n').split('\t')
        self.chrom = fields[0]
        self.source = fields[1]
        self.feature = fields[2]
        self.start = numstr(fields[3])
        self.end = numstr(fields[4])
        self.score = numstr(fields[5])
        self.strand = fields[6]
        self.frame = fields[7]
        self.attd = self.str_to_attd(fields[8])

    def str_to_attd(self, x):
        ret = OrderedDict()
        x = x.strip()
        if not x.endswith(';'): x += ';'
        for t in re.findall('(\S+)\s+"([\s\S]*?)";', x):
            ret[t[0]] = numstr(t[1])
        return ret

    def attd_to_str(self, x):
        """ Convert attribute dictionary to GTF string """
        return ' '.join(f'{k} "{v}";' for k, v in x.items())

    def __str__(self):
        return '\t'.join(map(str, [
            self.chrom, self.source, self.feature, self.start, self.end,
            self.score, self.strand, self.frame, self.attd_to_str(self.attd)
        ]))


class GFFRow(GTFRow):
    def str_to_attd(self, x):
        ret = OrderedDict()
        for kv in x.split(';'):
            k,v = kv.split('=')
            ret[k] = v
        return ret
    def attd_to_str(self, x):
        """ Convert attribute dictionary to GFF string """
        return ';'.join(f'{k}={v}' for k,v in x.items())

def numstr(x):
    """ Convert argument to numeric type.
    Attempts to convert argument to an integer. If this fails, attempts to
    convert to float. If both fail, return as string.
    Args:
        x: Argument to convert.
    Returns:
        The argument as int, float, or string.
    """
    try:
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return str(x)


def iterfile(filepath):
    """ Return filehandle for text or gzipped files"""
    with open(filepath, 'rb') as test_f:
        is_gzip = test_f.read(2) == b'\x1f\x8b'
    if is_gzip:
        return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'r')