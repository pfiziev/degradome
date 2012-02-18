from itertools import imap
import sys
from utils import *
import json
from pyfasta import *

hg19 = '/media/Data/UCLA/databases/genomes/hg19'


if __name__ == '__main__':

    if len(sys.argv) > 1:
        mapped_reads_015 = sys.argv[1]+ '.010anno.015map_hits'
        mapped_reads_020 = mapped_reads_015 + '.020stats'


    print 'input1:', mapped_reads_015
    print 'input2:', mapped_reads_020


    stats_type = 'bstats'
    pos_type = 'backward_position'

    # read the position statistics for stats_type
    stats = json.load(open(mapped_reads_020))[stats_type]


    # find out the position with the highest reads/regions ratio
    mpos = {}
    for r_type in stats:
        to_plot = sorted((int(k), float(stats[r_type][k][0])/stats[r_type][k][1]) for k in stats[r_type] if int(k) >= 0)[:200]
        mpos[r_type] = max(to_plot, key = lambda tp: tp[1])[0]

    RANGE_START = 0 # start the range relative to this position from mpos
    RANGE_END = 1   # end the range relative to this position from mpos
    mpos_ranges = dict((r_type, set(range(mpos[r_type] + RANGE_START, mpos[r_type] + RANGE_END))) for r_type in mpos)


    elapsed('reading position stats')

    # collect the regions that have reads in the mpos range
    regions = dict((r_type, []) for r_type in mpos)


    for reg in imap(json.loads, open(mapped_reads_015)):

        reg_peaks = mpos_ranges[reg['type']] & set(reg[pos_type])

        if len(reg_peaks) > 0:
            reg['peak'] = min(reg_peaks, key = lambda pos: abs(mpos[reg['type']] - pos))
            regions[reg['type']].append(reg)

    elapsed('reading position indices')


    # cut and output +-10 nucleotides from each region
    for r_type in regions:
        out = open(mapped_reads_020 + '.040_' + r_type + '.fa', 'w')
        chrom = None
        for reg in sorted(regions[r_type], key = lambda r: r['chrom']):
            if reg['chrom'] != chrom:
                chrom = reg['chrom']
                f = Fasta(os.path.join(hg19, chrom + '.fa'), flatten_inplace=True)

            cut_start = reg['end'] - reg['peak'] - 1
            cut_end = reg['end'] - reg['peak'] + 16

            out.write('>%s_%s_%d_%d_%s cut:%d-%d\n%s\n' % (r_type, chrom, reg['start'], reg['end'], reg['strand'], cut_start, cut_end,
                                                        f.sequence({'chr': chrom, 'start': cut_start, 'stop': cut_end, 'strand': reg['strand']})))



    json.dump(  {'regions' : regions, 'mpos' : mpos, 'RANGE_START' : RANGE_START, 'RANGE_END' : RANGE_END},
                open(mapped_reads_020+'.040regions', 'w'), indent = 1)
    elapsed('040get_regions_with_peaks')


