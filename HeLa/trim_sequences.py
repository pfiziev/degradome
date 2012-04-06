import sys


_rc_dict = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
def reverse_compl_seq(strseq):
    return ''.join(_rc_dict.get(c, c) for c in reversed(strseq))

if __name__ == '__main__':
    lno = 0
    entry = []
    discarded = 0
    total = 0
    for l in open(sys.argv[1]):
        l = l.strip()
        lno += 1
        entry.append(l)
        if lno == 4:
            total += 1
            if entry[0][0] != '@':
                print >> sys.stderr,  "ERROR:", ''.join(entry)
                exit(1)

            entry[1] = reverse_compl_seq(entry[1][:20])
            entry[3] = entry[3][:20][::-1]

            print '\n'.join(entry)
            
            entry = []
            lno = 0
 
