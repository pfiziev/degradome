import json
import os
from pprint import pformat
from utils import *

__author__ = 'pf'


if __name__ == '__main__':

    # read annotation
    anno = read_annotation()
    hit_reg = json.load(open(mapped_reads_010))
    elapsed('hit_reg annotation')
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


            region['forward_position'] = set()
            region['backward_position'] = set()
            region['chrom'] = chrom
            while tmp_ri < len(hits[chrom]) and region['end'] >= hits[chrom][tmp_ri][starti]:

                if hits[chrom][tmp_ri][strandi] == region['strand'] and \
                   region['type'] == hit_reg[hits[chrom][tmp_ri][sidi]] and \
                   partial_overlap(region['start'], region['end'], hits[chrom][tmp_ri][starti],  hits[chrom][tmp_ri][endi]):

                    region['forward_position'].add(       (hits[chrom][tmp_ri][starti] - region['start']) if region['strand'] == '+'
                                                     else (region['end'] - hits[chrom][tmp_ri][endi]))

                    region['backward_position'].add(      (region['end'] - hits[chrom][tmp_ri][starti]) if region['strand'] == '+'
                                                     else (hits[chrom][tmp_ri][endi] - region['start']))



#                    anno[chrom][tmp_ri]['overlap'] = 'full' if anno[chrom][tmp_ri]['start'] <= hit['start'] and hit['end'] <= anno[chrom][tmp_ri]['end'] else 'partial'

#                    hit_anno.append(anno[chrom][tmp_ri])
                tmp_ri += 1

            region['forward_position'] = sorted(region['forward_position'])
            region['backward_position'] = sorted(region['backward_position'])

            out.write(json.dumps(region)+'\n')
        elapsed(chrom)
    out.close()
    elapsed('annotate hits')







