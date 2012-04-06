__author__ = 'pf'
from utils import *
import sys

if __name__ == '__main__':
    # read annotation
#    anno = read_annotation()
#    anno_stats = {}
#    for chrom in anno:
#        anno_stats[chrom] = {}
#        for reg in anno[chrom]:
#            rt = reg['type']
#            if rt not in anno_stats[chrom]: anno_stats[chrom][rt] = 0
#            anno_stats[chrom][reg['type']] += 1
#
#    elapsed('anno_stats')

    if len(sys.argv) > 1:
        mapped_reads = sys.argv[1]
        mapped_reads_010 = sys.argv[1]+'.010anno'


    print 'mapped_reads:', mapped_reads
    print 'mapped_reads_010:', mapped_reads_010

    hit_reg = json.load(open(mapped_reads_010))
    elapsed('hit_reg annotation')
    read_chrom = {}
    for l in open(mapped_reads):
        sid, strand, chrom  = l.strip().split('\t')[:3]
        read_chrom[sid] = chrom
    elapsed('read_chrom')

    hit_stats = {}
    total = 0
    total_nonmit = 0
    for read, read_type in hit_reg.iteritems():
        rchrom = read_chrom[read]
        total += 1
        if rchrom == 'chrM':
            continue

        if read_type not in hit_stats: hit_stats[read_type] = 0
        hit_stats[read_type] += 1
        total_nonmit += 1
    elapsed('hit_stats')


    json.dump(hit_stats, open(mapped_reads+'.012read_stats', 'w'))


    print mapped_reads, total, total_nonmit

#
#hit_stats = {"exon_noncoding": 411441, "none": 269797, "exon": 3931354, "intron": 468056}
#from pylab import *
#
## make a square figure and axes
#figure(1, figsize=(6,6))
#ax = axes([0.1, 0.1, 0.8, 0.8])
#
#labels = 'Coding exons', 'Noncoding exons', 'Introns', 'Intergenic'
#fracs = [hit_stats[k] for k in ['exon', 'exon_noncoding', 'intron', 'none']]
#
#explode=(0, 0.05, 0, 0)
#pie(fracs, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True)
#title("HeLa", bbox={'facecolor':'0.8', 'pad':5})