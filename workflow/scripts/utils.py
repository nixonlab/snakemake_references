#! /usr/bin/env python
# -*- coding: utf-8 -*-

import re
from collections import OrderedDict


class GTFError(Exception):
    pass


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


def parse_gtf_attr(x):
    """ Convert attribute field to dictionary.
    Converts a GTF attribute string to a dictionary. The attribute string
    should be a semicolon-separated list of tag-value pairs. See
    http://www.ensembl.org/info/website/upload/gff.html.
    Args:
        x: Argument to convert.
    Returns:
        A dictionary containing the tag-value pairs.
    """
    if isinstance(x, dict):
        return x
    ret = OrderedDict()
    x = x.strip()
    if not x.endswith(';'): x += ';'
    for t in re.findall('(\S+)\s+"([\s\S]*?)";', x):
        ret[t[0]] = numstr(t[1])
    return ret


def attd_to_str(x):
    """ Convert dictionary to attribute string """
    return ' '.join('%s "%s";' % (k, v) for k,v in x.items())
