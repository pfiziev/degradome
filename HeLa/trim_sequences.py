import sys


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

            entry[1] = entry[1][:20]
            entry[3] = entry[3][:20]

            print '\n'.join(entry)
            
            entry = []
            lno = 0
 
