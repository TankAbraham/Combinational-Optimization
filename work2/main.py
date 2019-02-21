# -*- encoding:utf8 -*-

import time
import copy
import argparse
import numpy as np
from math import pow


def arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('n', help='value of n', type=int)
    parser.add_argument('b', help='value of b', type=float)
    args = parser.parse_args()
    return args.n, args.b


n, b = arg()


SAMPLE = 10
N = [i for i in range(n)]  # N
lenN = len(N)  # length of N
B = b  # b = 0.1
U = 100
ETA = 2
ENVIRONMENT = 'SLK'  # SLK or CON
# ENVIRONMENT = 'CON'  # SLK or CON
P={
    'a':[11, 13, 12, 15, 14],
    'v':[5, 3, 2, 9, 1],
    'omega':[4, 5, 2, 6, 3, 1]
}


# get_Omega() must be called after get_lambda()
def get_Omega():
    Omega = list()
    theLambda = P['lambda']
    for i in range(lenN):
        gross, thePower = 0, lenN - i - 2
        for j in range(lenN-1, i-1, -1):
            if thePower < 0: thePower = 0
            cur = theLambda[j]
            if j > 0 + i: cur = cur * B
            if j > 1 + i: cur = cur * pow((1+B), thePower)
            thePower = thePower - 1
            gross = gross + cur
        Omega.append(gross)
    return Omega


def get_median():
    omega = P['omega']
    if ENVIRONMENT is 'SLK':
        for i in range(len(omega)):
            if sum(omega[:i]) <= sum(omega[i:]):
                if sum(omega[:i+1]) >= sum(omega[i+1:]):
                    return i-1
    else:  # CON
        for i in range(len(omega)):
            if sum(omega[:i-1]) <= sum(omega[i-1:]):
                if sum(omega[:i]) >= sum(omega[i:]):
                    return i - 1


def get_lambda():
    k = P['median']
    omega = P['omega']
    l = list()
    if ENVIRONMENT is 'CON':
        for i in range(1, lenN+1):
            if i <= k:
                l.append(sum(omega[0:i]))
            else:
                l.append(sum(omega[i:]))
    else:  # SLK
        for i in range(1, lenN+1):
            if i <= k:
                l.append(sum(omega[0:i+1]))
            if i >= k + 1:
                if i <= lenN - 1:
                    l.append(sum(omega[i+1:]))
            if i == lenN:
                l.append(0)
    return l


def LB2(fixed, waiting):  # fixed, waiting are the 序号 of 工件
    X, Y, VY = list(), list(), list()
    for i in range(lenN):
        x, y = pow(P['Omega'][i], 1/(1+ETA)), pow(P['a'][i], ETA/(1+ETA))
        vy = y * P['v'][i]
        X.append(x), Y.append(y), VY.append(vy)
    part1 = 0
    for i, seq in enumerate(fixed):
        part1 = part1 + X[i] * Y[seq]
    for seq in sorted(fixed, reverse=True):
        Y.pop(seq)
    Y = sorted(Y)
    for j in range(len(Y)):
        part1 = part1 + X[j+i+1]*Y[j]
    part2 = 0
    for i, seq in enumerate(fixed):
        part2 = part2 + X[i] * VY[seq]
    for seq in sorted(fixed, reverse=True):
        VY.pop(seq)
    VY = sorted(VY)
    for j in range(len(VY)):
        part2 = part2 + X[j+i+1] * VY[j]
    return pow(U, -ETA) * part1 * pow(part2, ETA)


def obj(seq):
    X, Y = list(), list()
    for i in range(lenN):
        X.append(pow(P['Omega'][i], 1/(1+ETA)))
        Y.append(pow(P['a'][i], ETA/(1+ETA)))
    part1, part2 = 0, 0
    for i, iy in enumerate(seq):
        part1 = part1 + X[i] * Y[iy]
        part2 = part2 + X[i] * Y[iy] * P['v'][iy]
    return pow(U, -ETA) * part1 * pow(part2, ETA)


def explode(node):
    for i,e in enumerate(node[1]):
        temp = copy.deepcopy(node)
        temp[0].append(e)  # add an item into fixed
        temp[1].pop(i)  # remove an item from waiting
        temp[2] = LB2(temp[0], temp[1])  # Low Bound 2

        yield temp


def brunch_bound():
    theHB, buffer, best = HB2(), list(), None
    buffer.append(copy.deepcopy([[], [i for i in range(lenN)], 0]))  # [fixed, waiting, LB]
    while len(buffer):
        buffer = sorted(buffer, key=lambda e: e[2])  # 先按照LB进行排序,下面再按照展开的节点数进行排序
        buffer = sorted(buffer, key=lambda e: len(e[0]), reverse=True)  # 这样一样可以保证展开节点数多的、LB小的排在前面
        cur = buffer.pop(0)
        if len(cur[0]) is 0:  # node of root
            for node in explode(cur):
                if node[2] > theHB[1]: continue
                buffer.append(node)
        elif len(cur[0]) is lenN:  # node of leaf
            if cur[2] > theHB[1]: continue
            if best is None: best = cur
            elif cur[2] <= best[2]: best = cur
        else:  # node of intermedia
            for node in explode(cur):
                if node[2] > theHB[1]: continue
                buffer.append(node)
    return best


def HB2():  # vector
    X, Y = list(), list()                   # step 1
    for i in range(lenN):
        x = (i, pow(P['Omega'][i], 1/(1+ETA)))
        y = (i, pow(P['a'][i], ETA/(ETA+1)))
        X.append(x)
        Y.append(y)
    Y = sorted(Y, key=lambda e:e[1])
    seq = [(X[i][0], Y[i][0]) for i in range(lenN)]
    seq = sorted(seq, key=lambda e:e[0])  # sorted by the sequence of x
    seq1 = [e[1] for e in seq]   # sequence of step 1
    part1 = 0
    for i in range(lenN):
        part1 = part1 + X[i][1] * Y[i][1]
    part2 = 0
    for i,t in enumerate(Y):
        iY, valY = t
        part2 = part2 + X[i][1] * valY * P['v'][iY]
    obj1 = pow(U, -ETA) * part1 * pow(part2, ETA)

    Y = list()
    for i in range(lenN):                   # step 2
        Y.append((i, pow(P['a'][i], ETA / (ETA + 1)) * P['v'][i]))
    Y = sorted(Y, key=lambda e: e[1])  # VY
    seq = [(X[i][0], Y[i][0]) for i in range(lenN)]
    seq = sorted(seq, key=lambda e: e[0])  # sorted by the sequence of x
    seq2 = [e[1] for e in seq]  # sequence of step 2
    part2 = 0
    for i in range(lenN):
        part2 = part2 + X[i][1] * Y[i][1]
    Y = list()
    for i in range(lenN):                   # step 2
        y = (i, pow(P['a'][i], ETA / (ETA + 1)))
        Y.append(y)
    Y = sorted(Y, key=lambda e:e[0])
    part1 = 0
    for ix,iy in seq:
        part1 = part1 + X[ix][1] * Y[iy][1]
    obj2 = pow(U, -ETA) * part1 * pow(part2, ETA)

    if obj1 < obj2:                        # setp 3
        return (seq1, obj1)
    else:
        return (seq2, obj2)


if __name__ == '__main__':
    algo3, algo4, error_percentage = list(), list(), list()
    for i in range(SAMPLE):
        P['a'] = [np.random.randint(1, 100) for i in range(lenN)]
        P['alpha'] = [np.random.uniform(-0.5, -0.1) for i in range(lenN)]
        P['v'] = [np.random.randint(1, 100) for i in range(lenN)]
        P['omega'] = [np.random.randint(1, 100) for i in range(lenN + 1)]  # Omega and omega are different!
        P['median'] = get_median()
        P['lambda'] = get_lambda()
        P['Omega'] = get_Omega()  # Omega and omega are different!

        begin = time.time()
        Z_star = brunch_bound()[2]
        end = time.time()
        print("B&B running %.4f seconds" % (end-begin))
        algo4.append(end-begin)

        begin = time.time()
        Z = HB2()[1]
        end = time.time()
        print("heuristic running %.4f seconds" % (end-begin))
        algo3.append(end-begin)

        ep = (Z - Z_star) / Z_star
        print("error percentage is %f" % ep)
        error_percentage.append(ep)

    print("\nalgorithm 3 avg %f s" % (sum(algo3) / SAMPLE))
    print("algorithm 3 max %f s" % max(algo3))

    print("algorithm 4 avg %f s" % (sum(algo4) / SAMPLE))
    print("algorithm 4 max %f s" % max(algo4))

    print("error percentage avg %f s" % (sum(error_percentage) / SAMPLE))
    print("error percentage max %f s" % max(error_percentage))
