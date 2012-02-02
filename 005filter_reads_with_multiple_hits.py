from utils import *


if __name__ == '__main__':

    unique = set()
    for l in open(os.path.join('U87','U87.fq.noAdapters.contigFiltered.n0m1k1b.bt')):
        unique.add(l.split('\t')[0])
    elapsed('reading unique hits')

    out = open(os.path.join('U87','U87.fq.noAdapters.contigFiltered.non_unique.bt'),'w')
    for l in open(os.path.join('U87','U87.fq.noAdapters.contigFiltered.n0m20k20b.bt')):
        if l.split('\t')[0] not in unique:
            out.write(l)

    elapsed('done')
