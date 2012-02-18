import sys
from find_regions_with_peaks_shared import find_regions_with_peaks
from utils import mapped_reads_015, mapped_reads_020
from pfutils import *


snoRNA_annotation = '/home/pf/UCLA/databases/annotations/hg19_snoRNA'

if __name__ == '__main__':

    snoRNA = parsetsv(snoRNA_annotation,
                      ['_', 'chrom', ('start', int), ('end', int), 'name', '_', 'strand', '_', '_', 'type'],
                      groupby = 'chrom',
                      sortby =  'start')

    if len(sys.argv) > 1:
        mapped_reads_015 = sys.argv[1]+ '.010anno.015map_hits'

    mapped_reads_020 = mapped_reads_015 + '.020stats'
    mapped_reads_050 = mapped_reads_020 + '.050snoRNAs'

    regions = find_regions_with_peaks(mapped_reads_015, mapped_reads_020)
    elapsed('reading regions')
    for reg_type in regions:
        for reg in regions[reg_type]:
            reg['snoRNAs'] = [sr for sr in snoRNA.get(reg['chrom'],[]) if sr['strand'] == reg['strand'] and total_overlap(reg['start'], reg['end'], sr['start'], sr['end'])]

    res = dict((rt, [r for r in regions[rt] if len(r['snoRNAs']) > 0]) for rt in regions)
    for rt in res:
        print rt, ':', len(res[rt])
    json.dump(res, open(mapped_reads_050,'w'), indent = 1)
    elapsed('done')

