import datetime

__author__ = 'pf'



stime = datetime.datetime.now()
def elapsed(msg = None):
    print "[%s]" % msg if msg is not None else "+", "Elapsed: " , datetime.datetime.now() - stime

def partial_overlap(s1, e1, s2, e2):
    return s1 <= e2 and s2 <= e1 #and min(e1, e2) - max(s1, s2) >= 10

def total_overlap(s1, e1, s2, e2):
    return s1 <= s2 and e1 >= e2 or s2 <= s1 and e2 >= e1

