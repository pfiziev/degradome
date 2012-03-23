import json
import os
from pprint import pformat
import sys
from utils import *
from collections import defaultdict
__author__ = 'pf'


if __name__ == '__main__':

    # read annotation
    anno = read_annotation()

    if len(sys.argv) > 1:
        mapped_reads = sys.argv[1]
        mapped_reads_010 = sys.argv[1]+'.010anno'
        mapped_reads_015 = mapped_reads_010 + '.015map_hits'


    print 'mapped_reads:', mapped_reads
    print 'mapped_reads_010:', mapped_reads_010

    hit_reg = json.load(open(mapped_reads_010))
    elapsed('hit_reg annotation')
    # read bowtie hits
    hits = {}
    for l in open(mapped_reads):
        sid, strand, chrom, pos, seq = l.strip().split('\t')[:5]
        pos = int(pos)
        if chrom not in hits: hits[chrom] = []
        hits[chrom].append((sid,
                            strand,
                            pos,
                            pos + len(seq) - 1))

    sidi, strandi, starti, endi = 0, 1, 2, 3

    for chrom in hits:
        hits[chrom] = sorted(hits[chrom], key = lambda reg: reg[starti])
    elapsed('read hits')

    # annotate hits
    out = open(mapped_reads_015, 'w')
    for chrom in sorted(anno):
        if chrom not in hits:
            print 'skipping:', chrom
            continue

        print 'annotating hits on', chrom
        print 'hits:', len(hits[chrom]), 'annotated regions:', len(anno[chrom])

        r_index = 0
        for region in anno[chrom]:

            while r_index < len(hits[chrom]) and region['start'] > hits[chrom][r_index][endi]:
                r_index += 1
            tmp_ri = r_index


            region['forward_position'] = defaultdict(lambda : 0)
            region['backward_position'] = defaultdict(lambda : 0)

            region['chrom'] = chrom
            while tmp_ri < len(hits[chrom]) and region['end'] >= hits[chrom][tmp_ri][starti]:

                if hits[chrom][tmp_ri][strandi] == region['strand'] and \
                   region['type'] == hit_reg[hits[chrom][tmp_ri][sidi]] and\
                   partial_overlap(region['start'], region['end'], hits[chrom][tmp_ri][starti],  hits[chrom][tmp_ri][endi]):

                    region['forward_position'][  (hits[chrom][tmp_ri][starti] - region['start']) if region['strand'] == '+'
                                                 else (region['end'] - hits[chrom][tmp_ri][endi])] += 1

                    region['backward_position'][ (region['end'] - hits[chrom][tmp_ri][starti]) if region['strand'] == '+'
                                                 else (hits[chrom][tmp_ri][endi] - region['start'])] += 1


#                    anno[chrom][tmp_ri]['overlap'] = 'full' if anno[chrom][tmp_ri]['start'] <= hit['start'] and hit['end'] <= anno[chrom][tmp_ri]['end'] else 'partial'

#                    hit_anno.append(anno[chrom][tmp_ri])
                tmp_ri += 1

            MIN_READS_PER_PEAK = 3

            region['forward_position'] = sorted(pos for pos, count in region['forward_position'].iteritems() if count >= MIN_READS_PER_PEAK)
            region['backward_position'] = sorted(pos for pos, count in region['backward_position'].iteritems() if count >= MIN_READS_PER_PEAK)

            out.write(json.dumps(region)+'\n')
        elapsed(chrom)
    out.close()
    elapsed('annotate hits')







