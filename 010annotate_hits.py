import json
import os
from pprint import pformat
from utils import *

__author__ = 'pf'

#mapped_reads = os.path.join('U87','U87.fq.noAdapters.contigFiltered.n0m20k20b.bt')
mapped_reads = os.path.join('U87','U87.fq.noAdapters.contigFiltered.n0m1k1b.bt')

if __name__ == '__main__':

    # read annotation
    anno = {}
    for l in open('hg19.merged.to.ensg.all.tx.03.18.2011.txt.with.genetypes.final.txt'):
        _buf = l.strip().split('\t')
        strand = _buf[2]
        chrom = _buf[1]
        e_sts, e_ens = map(int, filter(None, _buf[8].split(','))), map(int, filter(None, _buf[9].split(',')))
        if chrom not in anno: anno[chrom] = []
        for i in xrange(len(e_sts)):
            anno[chrom].append({'start' : e_sts[i],
                                'end' : e_ens[i],
                                'type' : 'exon',
                                'strand' : strand,
                                'tid' : _buf[0],
                                'number' : i + 1 if strand == '+' else len(e_sts) - i})

            if i < len(e_sts) - 1:
                anno[chrom].append({'start' : e_ens[i] + 1,
                                    'end' : e_sts[i+1] - 1,
                                    'type' : 'intron',
                                    'strand' : strand,
                                    'tid' : _buf[0],
                                    'number' : i + 1 if strand == '+' else len(e_sts) - i - 1})

    for chrom in anno:
        anno[chrom] = sorted(anno[chrom], key = lambda reg: reg['start'])
    elapsed('read annotation')

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
    out = open(mapped_reads+'.010anno', 'w')
    for chrom in sorted(hits)   :
        print 'annotating hits on', chrom
        print len(hits[chrom]), len(anno[chrom])
#        print pformat(hits[chrom][:5])
#        print pformat(anno[chrom][:5])
        r_index = 0
        for  hit in hits[chrom]:
            hit_anno = []
#            while r_index < len(anno[chrom]) and not (anno[chrom][r_index]['start'] <= hit['end'] <= anno[chrom][r_index]['end']):

            while r_index < len(anno[chrom]) and hit['start'] > anno[chrom][r_index]['end']:
                r_index += 1
            tmp_ri = r_index

#            if hit_index < 5:
#                print '*' * 20
#                print pformat(hit)
#                print pformat(anno[chrom][r_index])
#                print r_index

            while tmp_ri < len(anno[chrom]) and hit['end'] >= anno[chrom][tmp_ri]['start']:

                if partial_overlap(hit['start'], hit['end'], anno[chrom][tmp_ri]['start'],  anno[chrom][tmp_ri]['end']):
                    anno[chrom][tmp_ri]['position'] = hit['start'] - anno[chrom][tmp_ri]['start'] if anno[chrom][tmp_ri]['strand'] == '+' \
                                                        else anno[chrom][tmp_ri]['end'] - hit['end']

                    anno[chrom][tmp_ri]['overlap'] = 'full' if hit['end'] <= anno[chrom][tmp_ri]['end'] else 'partial'

                    hit_anno.append(anno[chrom][tmp_ri])
                tmp_ri += 1

            hit['anno'] = hit_anno

            out.write(json.dumps(hit)+'\n')
        elapsed(chrom)
    out.close()
    elapsed('annotate hits')






        
