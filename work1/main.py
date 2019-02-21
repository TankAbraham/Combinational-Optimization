# -*- encoding:utf8 -*-

import time
import copy
from math import pow
import numpy as np
from scipy.optimize import linear_sum_assignment
import argparse


def arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('n', help='number of n', type=int)
    args = parser.parse_args()
    return args.n


SAMPLE = 10
N = [i for i in range(arg())]  # N
lenN = len(N)
M1 = np.zeros(shape=(lenN, lenN))
M2 = np.zeros(shape=(lenN, lenN))
overline_U = 200
ETA = 1.5
ENVIRONMENT = 'SLK'  # SLK or CON
# ENVIRONMENT = 'CON'  # SLK or CON
P={
    'p':[11, 13, 12, 15, 14],
    'alpha':[-0.32, -0.22, -0.35, -0.45, -0.25],
    'v':[5, 3, 2, 9, 1],
    'omega':[4, 5, 2, 6, 3, 1]
}


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


def assignment(fixed, waiting, m):
    m = np.delete(m, fixed, axis=0)  # axis = 0 is delete a specific row
    m = np.delete(m, [i for i in range(len(fixed))], axis=1)  # axis = 1 is delete a specific col
    row_idx, col_idx = linear_sum_assignment(m)
    cost = m[row_idx, col_idx].sum()
    return cost


def LB(fixed, waiting):
    l = list()
    for i,e in enumerate(fixed):
        lam = P['lambda'][i]
        ppi = P['p'][e]
        alpha = P['alpha'][e]
        f1 = pow(lam, 1/(1+ETA))
        f2 = pow(i+1, alpha)
        f3 = ppi * f2
        f4 = pow(f3, ETA/(ETA+1))
        l.append(f1 * f4)
    formula1 = sum(l)
    formula2 = assignment(fixed, waiting, M1)

    l = list()
    for i, e in enumerate(fixed):
        lam = P['lambda'][i]
        ppi = P['p'][e]
        alpha = P['alpha'][e]
        v = P['v'][e]  # e表示由工件所控制
        f1 = pow(lam, 1 / (1 + ETA))
        f2 = pow(i + 1, alpha)
        f3 = ppi * f2
        f4 = pow(f3, ETA / (ETA + 1))
        l.append(f1 * f4 * v)
    formula3 = sum(l)
    formula4 = assignment(fixed, waiting, M2)
    return pow(overline_U, -ETA) * (formula1 + formula2) * pow((formula3 + formula4), ETA)


def init_assignment():
    for r in range(lenN):
        for c in range(lenN):
            lam = P['lambda'][c]
            pi = P['p'][r]
            i = c + 1
            alpha = P['alpha'][r]
            v = P['v'][r]
            f1 = pow(lam, 1 / (1 + ETA))
            f2 = pi * pow(i, alpha)
            f3 = pow(f2, ETA / (ETA + 1))
            M1[r][c] = f1 * f3
            M2[r][c] = M1[r][c] * v


def explode(node, buffer):
    for i,e in enumerate(node[1]):
        temp = copy.deepcopy(node)
        temp[0].append(e)  # add an item into fixed
        temp[1].pop(i)  # remove an item from waiting
        temp[2] = LB(temp[0], temp[1])  # Low Bound 2
        buffer.append(temp)


def brunch_bound():
    theHB, buffer = HB()[2], list()
    best = [[], [], theHB]
    buffer.append(copy.deepcopy([[], [i for i in range(lenN)], 0]))  # [fixed, waiting, LB]
    while len(buffer):
        buffer = sorted(buffer, key=lambda e: e[2])  # 先按照LB进行排序,下面再按照展开的节点数进行排序
        buffer = sorted(buffer, key=lambda e: len(e[0]), reverse=True)  # 这样一样可以保证展开节点数多的、LB小的排在前面
        cur = buffer.pop(0)
        if len(cur[0]) is 0:  # node of root
            explode(cur, buffer)
        elif len(cur[0]) is lenN:  # node of leaf
            if cur[2] <= theHB and cur[2] <= best[2]:
                best = cur
        else:  # node of intermedia
            if cur[2] <= theHB:
                explode(cur, buffer)
    return best


def HB():
    row_idx, col_idx = linear_sum_assignment(M1)
    f1 = M1[row_idx, col_idx].sum()
    f2 = M2[row_idx, col_idx].sum()
    cost1 = pow(overline_U, -ETA) * f1 * pow(f2, ETA)

    row_idx, col_idx = linear_sum_assignment(M2)
    f1 = M1[row_idx, col_idx].sum()
    f2 = M2[row_idx, col_idx].sum()
    cost2 = pow(overline_U, -ETA) * f1 * pow(f2, ETA)

    print("high bound is %.4f and %.4f" % (cost1, cost2))
    return (cost1, cost2, min(cost1, cost2))


if __name__ == '__main__':
    algo3, algo4, error_percentage = list(), list(), list()
    for i in range(SAMPLE):
        P['p'] = [np.random.randint(1, 100) for i in range(lenN)]
        P['alpha'] = [np.random.uniform(-0.5, -0.1) for i in range(lenN)]
        P['v'] = [np.random.randint(1, 100) for i in range(lenN)]
        P['omega'] = [np.random.randint(1, 100) for i in range(lenN + 1)]
        P['median'] = get_median()
        P['lambda'] = get_lambda()
        init_assignment()

        begin = time.time()
        Z_star = brunch_bound()
        end = time.time()
        print("B&B running %.4f seconds" % (end-begin))
        algo4.append(end-begin)

        begin = time.time()
        Z = HB()[1]
        end = time.time()
        print("heuristic running %.4f seconds" % (end-begin))
        algo3.append(end-begin)

        ep = (Z - Z_star[2]) / Z_star[2]
        print("error percentage is %f" % ep)
        error_percentage.append(ep)

    print("\nalgorithm 3 avg %f s" % (sum(algo3) / SAMPLE))
    print("algorithm 3 max %f s" % max(algo3))

    print("algorithm 4 avg %f s" % (sum(algo4) / SAMPLE))
    print("algorithm 4 max %f s" % max(algo4))

    print("error percentage avg %f s" % (sum(error_percentage) / SAMPLE))
    print("error percentage max %f s" % max(error_percentage))
