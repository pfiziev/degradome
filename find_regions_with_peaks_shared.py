from itertools import imap
import json
from pfutils import *
from itertools import *


def find_regions_with_peaks(mapped_reads_015, mapped_reads_020):
    stats_type = 'bstats'
    pos_type = 'backward_position'

    # read the position statistics for stats_type
    stats = json.load(open(mapped_reads_020))[stats_type]


    # find out the position with the highest reads/regions ratio
    mpos = {}
    for r_type in stats:
        to_plot = sorted((int(k), float(stats[r_type][k][0])/stats[r_type][k][1]) for k in stats[r_type] if int(k) >= 0)[:200]
        mpos[r_type] = max(to_plot, key = lambda tp: tp[1])[0]

    RANGE_START = -3 # start the range relative to this position from mpos
    RANGE_END = 5   # end the range relative to this position from mpos
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

    return regions

