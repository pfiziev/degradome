from itertools import imap
import random
import sys
from find_regions_with_peaks_shared import find_regions_with_peaks
from utils import *
import json
from pyfasta import *

# this is made so it will work under linux and windows due to different path separators
#hg19 = reduce(os.path.join, '../../../databases/genomes/hg19'.split('/'))
hg19 = 'hg19'


if __name__ == '__main__':

    if len(sys.argv) > 1:
        mapped_reads_015 = sys.argv[1]+ '.010anno.015map_hits'#+'_min2'
        mapped_reads_020 = mapped_reads_015 + '.020stats'


    print 'input1:', mapped_reads_015
    print 'input2:', mapped_reads_020


    stats_type = 'bstats'
    pos_type = 'backward_position'

    regions = find_regions_with_peaks(mapped_reads_015, mapped_reads_020)

#    # read the position statistics for stats_type
#    stats = json.load(open(mapped_reads_020))[stats_type]
#
#
#    # find out the position with the highest reads/regions ratio
#    mpos = {}
#    for r_type in stats:
#        to_plot = sorted((int(k), float(stats[r_type][k][0])/stats[r_type][k][1]) for k in stats[r_type] if int(k) >= 0)[:200]
#        mpos[r_type] = max(to_plot, key = lambda tp: tp[1])[0]
#
#    RANGE_START = 0 # start the range relative to this position from mpos
#    RANGE_END = 1   # end the range relative to this position from mpos
#    mpos_ranges = dict((r_type, set(range(mpos[r_type] + RANGE_START, mpos[r_type] + RANGE_END))) for r_type in mpos)
#
#
#    elapsed('reading position stats')
#
#    # collect the regions that have reads in the mpos range
#    regions = dict((r_type, []) for r_type in mpos)
#
#
#    for reg in imap(json.loads, open(mapped_reads_015)):
#
#        reg_peaks = mpos_ranges[reg['type']] & set(reg[pos_type])
#
#        if len(reg_peaks) > 0:
#            reg['peak'] = min(reg_peaks, key = lambda pos: abs(mpos[reg['type']] - pos))
#            regions[reg['type']].append(reg)
#
#    elapsed('reading position indices')


    # cut and output +-10 nucleotides from each region
    seen = set()
    for r_type in regions:
        out = open(mapped_reads_020 + '.040_' + r_type + '.fa', 'w')
#        out_random = open(mapped_reads_020 + '.040_random_' + r_type + '.fa', 'w')
        chrom = None
        for i, reg in enumerate(sorted(regions[r_type], key = lambda r: r['chrom'])):
            if reg['chrom'] != chrom:
                chrom = reg['chrom']
                f = Fasta(os.path.join(hg19, chrom + '.fa'), flatten_inplace=True)

            if reg['strand'] == '+':
                cut_start = reg['end'] - reg['peak'][0] - 2
                cut_end = reg['end'] - reg['peak'][0] + 30
            else:
                cut_end = reg['start'] + reg['peak'][0] + 2
                cut_start = reg['start'] + reg['peak'][0] - 30


#            random_peak = random.randint(0, reg['end'] - reg['start'])
#            random_cut_start = reg['end'] - random_peak - 2
#            random_cut_end   = reg['end'] - random_peak + 30

            if hash_region(reg) in seen:
                continue
            seen.add(hash_region(reg))


#            out.write('>%s_%d %s\n%s\n' % (r_type, i, json.dumps({
#                                                                'peak' : reg['peak'],
#                                                                'region' : '%(chrom)s:%(start)d-%(end)d' % reg,
#                                                                'strand': reg['strand'],
#                                                                'tid' : reg.get('tid'),
#                                                                'type' : reg['type']
#                                                                }),
#                                                    f.sequence({'chr': chrom,
#                                                                'start': reg['start'],
#                                                                'stop': reg['end'],
#                                                                'strand': reg['strand']})))

            out.write('>%s_%d %s\n%s\n' % (r_type, i, json.dumps({
                                                                'peak' : reg['peak'],
                                                                'region' : '%(chrom)s:%(start)d-%(end)d' % reg,
                                                                'strand': reg['strand'],
                                                                'tid' : reg.get('tid'),
                                                                'type' : reg['type'],
                                                                'cut': [cut_start, cut_end]}),
                                                    f.sequence({'chr': chrom, 'start': cut_start, 'stop': cut_end, 'strand': reg['strand']})))


#            out_random.write('>%s_%d %s\n%s\n' % (r_type, i, json.dumps({
#                                                                'peak' : random_peak,
#                                                                'region' : '%(chrom)s:%(start)d-%(end)d' % reg,
#                                                                'strand': reg['strand'],
#                                                                'tid' : reg.get('tid'),
#                                                                'type' : reg['type'],
#                                                                'cut': [random_cut_start, random_cut_end]}),
#                                                    f.sequence({'chr': chrom, 'start': random_cut_start, 'stop': random_cut_end, 'strand': reg['strand']})))



#    json.dump(  {'regions' : regions, 'mpos' : mpos, 'RANGE_START' : RANGE_START, 'RANGE_END' : RANGE_END},
#                open(mapped_reads_020+'.040regions', 'w'), indent = 1)
    json.dump(regions, open(mapped_reads_020+'.040regions', 'w'), indent = 1)
    elapsed('040get_regions_with_peaks')

## COMPUTE OVERLAP BETWEEN DIFFERENT DATASETS
#
#
#import json
#hela = json.load(open('/mnt/hoffman2/u/home/mcdb/pfiziev/projects/degradome/HeLa/ALL.fastq.trimmed.n0m1k1b.bt.010anno.015map_hits.020stats.040regions'))
#u87 = json.load(open('/mnt/hoffman2/u/home/mcdb/pfiziev/projects/degradome/U87/U87.fq.noAdapters.n0m1k1b.bt.010anno.015map_hits.020stats.040regions'))
#u87_dkd = json.load(open('/mnt/hoffman2/u/home/mcdb/pfiziev/projects/degradome/U87_DroshaKD/Dro.fq.003noAdapters.n0m1k1b.bt.010anno.015map_hits_min2.020stats.040regions'))
#def rkey(reg):
#    return "%(chrom)s:%(start)d-%(end)d" % reg
#
#rtype = 'exon_noncoding'
#venn([map(rkey, hela[rtype]), map(rkey, u87[rtype]), map(rkey, u87_dkd[rtype])], ["HeLa", "U87", "U87_DroshaKD"], fill="number", show_names=True)
#
