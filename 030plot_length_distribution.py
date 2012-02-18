from collections import defaultdict
import json
import random
import sys
from utils import *

__author__ = 'pf'

#input = 'U87/x'

if __name__ == '__main__':

    # take into account hits with start positions between hit_start and hit_end
    hit_start = 20
    hit_end = 50
    reg_type = 'intron'
    pos_type = 'backward_position'


    seen = set()

    if len(sys.argv) > 1:
        mapped_reads_015 = sys.argv[1]+ '.010anno' + '.015map_hits'

    target_regions = []
    all_regions = []

    print 'input:', mapped_reads_015

    for l in open(mapped_reads_015):
        reg = json.loads(l)
        key = '%(chrom)s %(start)d %(end)d' % reg
        if reg['type'] != reg_type or key in seen: continue
        seen.add(key)

        reg_length = reg['end'] - reg['start'] + 1

        if reg['type'] == 'exon' and reg_length > 20000:
            continue

        if any(hit_start <= pos <= hit_end for pos in reg[pos_type]):
            target_regions.append(reg_length)

        all_regions.append(reg_length)

    random_regions = random.sample(all_regions, 1000)


    elapsed('initial stats')



    # plot the histogram


    import matplotlib.pyplot as plt


    # first we'll do it the default way, with gaps on weekends
    fig = plt.figure()
    ax = fig.add_subplot(111)

    #    ax.plot([c[0] for c in to_plot], [c[1] for c in to_plot], 'o')
    ax.hist(target_regions, 100, facecolor='green', alpha=0.75)
    ax.hist(random_regions, 100, facecolor='blue', alpha=0.75)
    fig.savefig('030length_histogram.png')

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
