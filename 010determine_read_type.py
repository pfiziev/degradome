import json
import os
from pprint import pformat
import sys
from utils import *

__author__ = 'pf'


if __name__ == '__main__':

    # read annotation
    anno = read_annotation()

    # read bowtie hits
    hits = {}
    if len(sys.argv) > 1:
        mapped_reads = sys.argv[1]
        mapped_reads_010 = mapped_reads + '.010anno'

    print 'input:', mapped_reads

    for l in open(mapped_reads):
        try:
            sid, strand, chrom, pos, seq = l.strip().split('\t')[:5]
        except:
            print '|'+l.strip()+'|'
            raise
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






        
