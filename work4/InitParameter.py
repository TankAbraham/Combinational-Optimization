# encoding:utf8

import sys
from random import randint


def parameter(n, f):
    while True:
        l = [randint(1, int(2*n/f)+1) for i in range(f)]
        if sum(l) is n:
            break
    g = list()
    for i in range(f):
        d = {
            'idx': i,
            'p': [randint(1, 100) for j in range(l[i])],
            'xi': [randint(1, 100) for j in range(l[i])],
            'omega': randint(1, 100) ,
            's': randint(1, 100),
            'n': l[i]
        }
        g.append(d)
    return g


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: python3 InitParameter.py n f instance')
        exit(1)
    n, f, instance = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    fname = 'data/%d_%d.txt' % (n, f)
    fout = open(fname, mode='w', encoding='utf8')
    for i in range(instance):
        fout.write("%s\n" % str(parameter(n, f)))
