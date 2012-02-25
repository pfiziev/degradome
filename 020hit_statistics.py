from collections import defaultdict
import json
import sys
from utils import *

__author__ = 'pf'

#input = 'U87/x'

if __name__ == '__main__':


    seen = set()
    fstats = defaultdict(lambda: defaultdict(lambda: 0))
    bstats = defaultdict(lambda: defaultdict(lambda: 0))
    lstats = defaultdict(lambda: defaultdict(lambda: 0))

    if len(sys.argv) > 1:
        mapped_reads_015 = sys.argv[1]+ '.010anno' + '.015map_hits'

    print 'input:', mapped_reads_015

    for l in open(mapped_reads_015):
        reg = json.loads(l)
        key = "%(chrom)s %(start)d %(end)d" % reg
        if key in seen: continue
        seen.add(key)

        if reg['type'] == 'exon' and reg['end'] - reg['start'] + 1 > 20000:
            continue

        for pos in reg['forward_position']:
            fstats[reg['type']][pos] += 1

        for pos in reg['backward_position']:
            bstats[reg['type']][pos] += 1

        # count the regions that have (end - start) valid positions i.e. have length of (end - start + 1).
        # remember, end is the 0-based position of the last nucleotide in the region ( seq[end] is a valid nucleotide position)
        lstats[reg['type']][reg['end'] - reg['start']] += 1

    def amend_stats(stats):
        for regtype in stats:
            regs_lengths = sorted(lstats[regtype], reverse = True)
            regs_index = 0
            regs_greater = 0
            for pos in sorted(stats[regtype], reverse = True):

                while regs_index < len(regs_lengths) and regs_lengths[regs_index] >= pos:
                    regs_greater += lstats[regtype][regs_lengths[regs_index]]
                    regs_index += 1

                stats[regtype][pos] = [stats[regtype][pos], regs_greater]



    elapsed('initial stats')

    amend_stats(fstats)
    elapsed('ammend_stats: fstats')

    amend_stats(bstats)
    elapsed('ammend_stats: bstats')


    json.dump({'fstats': fstats, 'bstats' : bstats}, open(mapped_reads_015+'.020stats','w'), indent = 1)
    elapsed('020hit_statistics.py')

    """
    # to plot the data

import matplotlib.pyplot as plt
import json

stats = json.load(open('/home/pf/UCLA/WINTER 2012/Grace/degradome/U87_DroshaKD/Dro.fq.003noAdapters.n0m1k1b.bt.010anno.015map_hits.020stats'))
rtype = 'intron'
fstats = stats['bstats']
to_plot = sorted((int(k)+1, float(fstats[rtype][k][0])/fstats[rtype][k][1]) for k in fstats[rtype] if int(k) >= 0)[:200]


# first we'll do it the default way, with gaps on weekends
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot([c[0] for c in to_plot], [c[1] for c in to_plot], 'o')

    """

#    # read annotation
#    anno = read_annotation()
#
#    reg_stats = defaultdict(lambda: defaultdict(lambda: 0))
#
#    for chrom in anno:
#        for reg in anno[chrom]:
#            reg_len = reg['end'] - reg['start']
#
#            if reg_len < 30: continue
#
#            reg_stats[reg['type']][reg_len ] += 1
#
#
##    logf(reg_stats)
#    elapsed('region stats')
#
#    # reads stats
#    stats = defaultdict(lambda: defaultdict(lambda: 0))
#    reg_with_hits_stats = defaultdict(lambda: defaultdict(lambda: 0))
#    _seen_regs = set()
#    for l in open(input):
#        mr = json.loads(l)
#        rtype = 'exon' if any(a['type'] == 'exon' for a in mr['anno']) else 'intron'
#
#        for anno in mr['anno']:
#            if anno['strand'] != mr['strand'] or anno['overlap'] != 'full' or anno['type'] != rtype: continue
#
#
#            stats[anno['type']][anno['position']] += 1
#
#            k = "%s %d %d" % (anno['tid'], anno['start'], anno['end'])
#            if k not in _seen_regs:
#                reg_with_hits_stats[anno['type']][anno['end'] - anno['start'] ] += 1
#                _seen_regs.add(k)
#
#    elapsed('hit stats')
#    # normalize stats
#    final_stats = defaultdict(list)
#
#    for reg in stats:
#        regs_greater = 0
#        regs_index = 0
#        regs_lengths = sorted(reg_stats[reg], reverse=True)
#
#        regs_wh_greater = 0
#        regs_wh_index = 0
#        regs_wh_lengths = sorted(reg_with_hits_stats[reg], reverse=True)
#
#
#
#
#        for pos in sorted(stats[reg], reverse = True):
#
#            while regs_index < len(regs_lengths) and regs_lengths[regs_index] > pos:
#                regs_greater += reg_stats[reg][regs_lengths[regs_index]]
#                regs_index += 1
#
#
#            while regs_wh_index < len(regs_wh_lengths) and regs_wh_lengths[regs_wh_index] > pos:
#                regs_wh_greater += reg_with_hits_stats[reg][regs_wh_lengths[regs_wh_index]]
#                regs_wh_index += 1
#
#
#
#
##            stats[reg][pos] = (stats[reg][pos], regs_greater)
#            final_stats[reg].append([pos, stats[reg][pos], regs_greater, regs_wh_greater])
##            stats[reg][pos] /= float(regs_greater)
#
#    open(input+'.020stats','w').write(json.dumps(final_stats))
#
#
#
