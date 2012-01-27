import json
import os
from pprint import pformat
from utils import *

__author__ = 'pf'


if __name__ == '__main__':

    # read annotation
    anno = read_annotation()

    # read bowtie hits
    hits = {}
    for l in open(mapped_reads):
        sid, strand, chrom, pos, seq, _, _ = l.strip().split('\t')
        pos = int(pos)
        if chrom not in hits: hits[chrom] = []
        hits[chrom].append((sid,
                            strand,
                            pos,
                            pos + len(seq)))

    for chrom in hits:
        hits[chrom] = sorted(hits[chrom], key = lambda reg: reg[2])
    elapsed('read hits')


    # determine the type of the region each read maps to
    hit_reg = {}
    out = open(mapped_reads_010, 'w')
    for chrom in sorted(hits)   :
        print 'annotating hits on', chrom
        print len(hits[chrom]), len(anno[chrom])

        r_index = 0
        for sid, strand, start, end in hits[chrom]:
            hit_anno = []

            while r_index < len(anno[chrom]) and start > anno[chrom][r_index]['end']:
                r_index += 1
            tmp_ri = r_index

            read_type = 'none'
            while tmp_ri < len(anno[chrom]) and end >= anno[chrom][tmp_ri]['start']:

                if anno[chrom][tmp_ri]['strand'] == strand and partial_overlap(start, end, anno[chrom][tmp_ri]['start'],  anno[chrom][tmp_ri]['end']):
                    if anno[chrom][tmp_ri]['type'] == 'exon':
                        read_type = 'exon'
                        break
                    else:
                        read_type = 'intron'

                tmp_ri += 1

            hit_reg[sid] = read_type


        elapsed(chrom)
    json.dump(hit_reg, out)
    out.close()
    elapsed('annotate hits')






        
