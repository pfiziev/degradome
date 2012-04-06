__author__ = 'pf'

# THIS WAS RUN IN DREAMPIE

#f = open('/mnt/hoffman2/u/home/mcdb/pfiziev/projects/degradome/peakData.10.hg19.U87.both.unique.lax')
#peaks = {}
#for l in f:
#    chrom, strand, start, end = l.strip().split(':')
#    strand = '+' if strand == '1' else '-'
#    if chrom not in peaks: peaks[chrom] = []
#    peaks[chrom].append((strand, int(start), int(end)))
#f.close()
#
#def partial_overlap(s1, e1, s2, e2):
#    return s1 <= e2 and s2 <= e1 #and min(e1, e2) - max(s1, s2) >= 10
#
#stats = {}
#import json
#regs = json.load(open('/mnt/hoffman2/u/home/mcdb/pfiziev/projects/degradome/U87/U87.fq.noAdapters.n0m1k1b.bt.010anno.015map_hits.020stats.040regions'))
#for reg_type in regs:
#    stats[reg_type] = 0
#    for reg in regs[reg_type]:
#        for p_str, p_st, p_en in peaks.get(reg['chrom'],[]):
#            if p_str == reg['strand'] and partial_overlap(p_st, p_en, reg['start'], reg['end']):
#                stats[reg_type] += 1
#print stats