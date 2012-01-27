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
        hits[chrom].append({'sid' : sid,
                            'strand' : strand,
                            'start' : pos,
                            'end' : pos + len(seq),
                            'chrom' : chrom })

    for chrom in hits:
        hits[chrom] = sorted(hits[chrom], key = lambda reg: reg['start'])
    elapsed('read hits')

    # annotate hits
    out = open(mapped_reads_015, 'w')
    for chrom in sorted(anno)   :
        if chrom not in hits:
            print 'skipping:', chrom
            continue

        print 'annotating hits on', chrom
        print 'hits:', len(hits[chrom]), 'annotated regions:',len(anno[chrom])
        #        print pformat(hits[chrom][:5])
        #        print pformat(anno[chrom][:5])
        r_index = 0
        for region in anno[chrom]:
            #            while r_index < len(anno[chrom]) and not (anno[chrom][r_index]['start'] <= hit['end'] <= anno[chrom][r_index]['end']):

            while r_index < len(hits[chrom]) and region['start'] > hits[chrom][r_index]['end']:
                r_index += 1
            tmp_ri = r_index

            #            if hit_index < 5:
            #                print '*' * 20
            #                print pformat(hit)
            #                print pformat(anno[chrom][r_index])
            #                print r_index

            region['forward_position'] = set()
            region['backward_position'] = set()
            region['chrom'] = chrom
            while tmp_ri < len(hits[chrom]) and region['end'] >= hits[chrom][tmp_ri]['start']:

                if hits[chrom][tmp_ri]['strand'] == region['strand'] and \
                   region['type'] == hit_reg[hits[chrom][tmp_ri]['sid']] and \
                   partial_overlap(region['start'], region['end'], hits[chrom][tmp_ri]['start'],  hits[chrom][tmp_ri]['end']):

                    region['forward_position'].add(       (hits[chrom][tmp_ri]['start'] - region['start']) if region['strand'] == '+'
                                                     else (region['end'] - hits[chrom][tmp_ri]['end']))

                    region['backward_position'].add(      (region['end'] - hits[chrom][tmp_ri]['start']) if region['strand'] == '+'
                                                     else (hits[chrom][tmp_ri]['end'] - region['start']))



#                    anno[chrom][tmp_ri]['overlap'] = 'full' if anno[chrom][tmp_ri]['start'] <= hit['start'] and hit['end'] <= anno[chrom][tmp_ri]['end'] else 'partial'

#                    hit_anno.append(anno[chrom][tmp_ri])
                tmp_ri += 1

            region['forward_position'] = sorted(region['forward_position'])
            region['backward_position'] = sorted(region['backward_position'])

            out.write(json.dumps(region)+'\n')
        elapsed(chrom)
    out.close()
    elapsed('annotate hits')







