import datetime
import json
import os

__author__ = 'pf'

def strongtype(method):
    """
    Decorator for strong-typing methods.

    Usage:

    >>> @strongtype
    ... def foo(a, b = int, c = float):
    ...     print a, type(a)
    ...     print b, type(b)
    ...     print c, type(c)
    ...
    ... foo(a ='1',  c= 3, b = 32)
    1 <type 'str'>
    32 <type 'int'>
    3.0 <type 'float'>

    """

    import inspect
    args, varargs, keywords, defaults = inspect.getargspec(method)

    nargs = len(args) - len(defaults) # number of normal arguments in the method definition

    types = dict((_kw, _type) for _kw, _type in zip(args[nargs:], defaults)) # a dictionary holding the (arg_name, type) mapping

    def new_method(*cur_args, **cur_kwargs):

        # assume that if cur_args contains more elements than nargs then the remainder of it (cur_args[nargs:]) contains keyword arguments
        # that were not stated explicitly. For example, if the definition of the method is "def foo(a, b = int)" a call like "foo('1','2')" would map the
        # value '2' to the definition of argument "b".
        cur_args = list(cur_args[:nargs]) + [types[_kw](_val) for _kw, _val in zip(args[nargs:], cur_args[nargs:])]

        method(*cur_args,
            **dict((_kw, types.get(_kw, lambda x: x)(_val)) # cast the keyword arguments
            for _kw, _val in cur_kwargs.items()))

    return new_method

#mapped_reads = os.path.join('U87','U87.fq.noAdapters.contigFiltered.n0m20k20b.bt')
mapped_reads = os.path.join('U87','U87.fq.noAdapters.contigFiltered.n0m1k1b.bt')
mapped_reads_010 = mapped_reads + '.010anno'
mapped_reads_015 = mapped_reads + '.015map_hits'

logfile = open(os.path.join('U87','log'), 'w')
def logf(msg):
    logfile.write((msg if isinstance(msg, str) else json.dumps(msg)) +'\n')

global_stime = datetime.datetime.now()
def elapsed(msg = None):
    print "[%s]" % msg if msg is not None else "+", "Last:" , datetime.datetime.now() - elapsed.stime, '\tTotal:', datetime.datetime.now() - global_stime
    elapsed.stime = datetime.datetime.now()

elapsed.stime = datetime.datetime.now()


def partial_overlap(s1, e1, s2, e2):
    return s1 <= e2 and s2 <= e1 #and min(e1, e2) - max(s1, s2) >= 10

def total_overlap(s1, e1, s2, e2):
    return s1 <= s2 and e1 >= e2 or s2 <= s1 and e2 >= e1


def read_annotation():
    anno = {}
    min_reg_length = 25

    for l in open('hg19.merged.to.ensg.all.tx.03.18.2011.txt.with.genetypes.final.txt'):
        _buf = l.strip().split('\t')
        strand = _buf[2]
        chrom = _buf[1]
        e_sts, e_ens = map(int, filter(None, _buf[8].split(','))), map(int, filter(None, _buf[9].split(',')))
        if chrom not in anno: anno[chrom] = []
        for i in xrange(len(e_sts)):
            if e_ens[i] - e_sts[i] >= min_reg_length:
                anno[chrom].append({'start' : e_sts[i],
                                    'end' : e_ens[i],
                                    'type' : 'exon',
                                    'strand' : strand,
                                    'tid' : _buf[0],
                                    'number' : i + 1 if strand == '+' else len(e_sts) - i})

            if i < len(e_sts) - 1 and e_sts[i+1] - e_ens[i] - 2 >= min_reg_length:
                anno[chrom].append({'start' : e_ens[i] + 1,
                                    'end' : e_sts[i+1] - 1,
                                    'type' : 'intron',
                                    'strand' : strand,
                                    'tid' : _buf[0],
                                    'number' : i + 1 if strand == '+' else len(e_sts) - i - 1})

    for chrom in anno:
        anno[chrom] = sorted(anno[chrom], key = lambda reg: reg['start'])

    elapsed('read annotation')

    return anno