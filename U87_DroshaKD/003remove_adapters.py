import sys
from utils import *


DRO_PRIMER = 'AGCGCTAGTC'
DRO_INPUT = os.path.join('U87_DroshaKD','Dro.fq')
DRO_OUT = open(DRO_INPUT+'.003noAdapters','w')

if __name__ == '__main__':
    lno = 0
    entry = []
    discarded = 0
    total = 0
    for l in open(DRO_INPUT):
        lno += 1
        entry.append(l)
        if lno == 4:
            total += 1
            if entry[0][0] != '@':
                print >> sys.stderr,  "ERROR:", ''.join(entry)
                exit(1)

            parts = entry[1].split(DRO_PRIMER)

            if len(parts) == 2 and len(parts[0]) >= 20:
                entry[1] = parts[0]+'\n'
                entry[3] = entry[3][0:len(entry[1]) - 1]
                DRO_OUT.write(''.join(entry)+'\n')
            else:
                discarded += 1

            entry = []
            lno = 0
    DRO_OUT.close()
    elapsed('Total:' + str(total) + '\tDiscarded:'+str(discarded)+'\t'+str(float(100*discarded)/total))