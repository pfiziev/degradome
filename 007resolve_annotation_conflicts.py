from utils import *

INTERGENIC = 0
EXON = 1
INTRON = 2

rtypes = ['intergenic', 'exon', 'intron']

if __name__ == '__main__':
    anno = read_annotation_original()
    maxl = 0
    for chrom in anno:
        _maxl = max(r['end'] for r in anno[chrom])
        if _maxl > maxl: maxl = _maxl

        anno[chrom] = {'+' : sorted((r for r in anno[chrom] if r['strand'] == '+'), key = lambda reg: reg['start']),
                       '-' : sorted((r for r in anno[chrom] if r['strand'] == '-'), key = lambda reg: reg['start'])}

    bitmap = [INTERGENIC for i in xrange(maxl+1)]
    elapsed('reading annotation')
    print 'bitmap size:', len(bitmap)

    new_anno = {}
    seen = set()

    for chrom in anno:

        new_anno[chrom] = []

        for strand in anno[chrom]:

            # reset the bitmap
            for i in xrange(len(bitmap)) :
                bitmap[i] = INTERGENIC



            for reg in anno[chrom][strand]:

                # skip identical regions
                rkey = '%s %d %d' % (chrom, reg['start'], reg['end'])
                if rkey in seen: continue
                seen.add(rkey)


                if reg['type']== 'exon':
                    for i in xrange(reg['start'], reg['end']+1):
                        bitmap[i] = EXON

                if reg['type'] == 'intron':
                    for i in xrange(reg['start'], reg['end']+1):
                        if bitmap[i] == INTERGENIC: bitmap[i] = INTRON

            pos = 0
            reg_start, reg_end, reg_type = 0, 0, INTERGENIC
            while pos < len(bitmap):
                reg_type = bitmap[pos]
                reg_start = pos

                while pos < len(bitmap) and reg_type == bitmap[pos]:
                    pos += 1
                reg_end = pos - 1

                if reg_type != INTERGENIC:
                    new_anno[chrom].append({'type' : rtypes[reg_type], 'start': reg_start, 'end' : reg_end, 'strand' : strand, 'chrom' : chrom})

        elapsed(chrom)

    json.dump(new_anno, open('hg19.merged.to.ensg.all.tx.03.18.2011.txt.with.genetypes.final.txt.007noConflicts', 'w'), indent = 1)
    elapsed('dump')