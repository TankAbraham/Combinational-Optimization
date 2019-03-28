# encoding:utf8

import sys
from time import clock
from copy import deepcopy
from math import factorial


eta, alpha, beta, gamma, a, b = 2, 1, 1, 1, 0.1, 0.2
g = [
    {
        'idx' : 0,
        'p': [2, 4, 6],
        'xi': [2, 3, 3],
        'omega': 1,
        's': 2,
        'n': 3
    },
    # {  # temp for debug
    #     'idx' : 1,
    #     'p': [5, 8, 4],
    #     'xi': [3, 2, 2],
    #     'omega': 1,
    #     's': 4,
    #     'n': 3
    # },
    {  # real
        'idx' : 1,
        'p': [5, 8, 4, 3],
        'xi': [3, 2, 2, 3],
        'omega': 1,
        's': 4,
        'n': 4
    },
    {
        'idx' : 2,
        'p': [9, 4],
        'xi': [4, 2],
        'omega': 1,
        's': 3,
        'n': 2
    },
    {
        'idx' : 3,
        'p': [6, 3, 3, 4],
        'xi': [1, 1, 3, 2],
        'omega': 1,
        's': 5,
        'n': 4
    },
    {
        'idx': 4,
        'p': [8, 10, 6],
        'xi': [4, 3, 1],
        'omega': 1,
        's': 6,
        'n': 3
    }
]


def Z(f, gs):  # group sequence
    constant = eta**(-eta/(eta+1)) + eta**(1/(eta+1))
    p1 = 0
    for i in range(f):
        for j in range(gs[i]['n']):
            power = sum([gs[l]['n'] for l in range(i, f)])-j-1
            p11 = ((1+b)**(f-i-1) * (1+a)**power)**(1/(eta+1))
            p12 = gs[i]['pi'][j]**(eta/(eta+1))
            p1 = p1 + p11 * p12
    p2 = 0
    for i in range(f):
        power = sum([gs[l]['n'] for l in range(i, f)])
        p21 = ((1+b)**(f-i-1) * (1+a)**power)**(1/(eta+1))
        p22 = (gs[i]['omega'] * gs[i]['s'])**(eta/(eta+1))
        p2 = p2 + p21 * p22
    z = constant * (p1 + p2)
    return z


def derive(node):
    l = list()
    idx = node['i']
    times = factorial(len(node['s'][idx]))
    temp = deepcopy(node)
    for i in range(times):
        temp = deepcopy(temp)
        temp['i'] = node['i'] + 1
        i = i % len(node['s'][idx])
        try:
            temp['s'][idx][i], temp['s'][idx][i+1] = temp['s'][idx][i+1], temp['s'][idx][i]
        except IndexError:
            temp['s'][idx][i], temp['s'][idx][0] = temp['s'][idx][0], temp['s'][idx][i]
        l.append(temp)
    return l


def search(f, seq):  # seq was sorted by n
    s, node, last = list(), [seq[0]], seq[0]['n']
    for i in range(1, f):
        if seq[i]['n'] is last:
            node.append(seq[i])
        else:
            s.append(node)
            node = list()
            node.append(seq[i])
            last = seq[i]['n']
    s.append(node)
    buffer = [{'s':s, 'i':0}]
    best = {'z': Z(f, seq), 's':seq}
    while(len(buffer)):
        cur = buffer.pop(0)
        if cur['i'] is len(cur):  # leaf
            l = list()
            for i in range(len(cur['s'])):
                for j in range(len(cur['s'][i])):
                    l.append(cur['s'][i][j])
            z = Z(f, l)
            if z < best['z']:
                best['z'], best['s'] = z, l
            continue
        for s in derive(cur):
            buffer.append(s)
    return best['z'], best['s']


def heuristic(f):
    for i in range(f):  # step 1
        pi = [g[i]['xi'][j] * g[i]['p'][j] for j in range(g[i]['n'])]
        g[i]['pi'] = sorted(pi)
    g1 = sorted(g, key=lambda x:x['omega'] * x['s'])  # step 2
    g2 = sorted(g, key=lambda x:x['n'], reverse=True)  # step 3
    z1 = Z(f, g1)
    z2, g2 = search(f, g2)
    gs = deepcopy(g1) if z1 < z2 else deepcopy(g2)  # step 4
    z = deepcopy(z1) if z1 < z2 else deepcopy(z2)
    target = {'g' : gs, 'z' : z}
    for i in range(f-1):  # step 5,6,7,8,9
        cur = deepcopy(target)
        for j in range(i+1, f):
            temp = deepcopy(cur)
            temp['g'][i], temp['g'][j] = temp['g'][j], temp['g'][i]
            temp['z'] = Z(f, temp['g'])
            if temp['z'] < target['z']:
                target = deepcopy(temp)
    return target


# set.      p = p1 + p2
#           p1 = p11 + p12, p2 = p21 + p22
#           p12 = p121 + p122, p22 = p221 + p222
# hence.    p = p11(p121+p122)+p21(p221+p222)
def LB(fix, wait):
    wait = sorted(wait, key=lambda x:x['omega']*x['s'])
    c = sorted([wait[i]['n'] for i in range(len(wait))], reverse=True)
    p11 = (eta**(-eta/(eta+1)) + eta**(1/(eta+1))) * alpha**(1/(eta+1)) * beta**(eta/(eta+1))
    p21 = (eta**(-eta/(eta+1)) + eta**(1/(eta+1))) * alpha**(1/(eta+1)) * gamma**(eta/(eta+1))
    p121 = 0
    for i in range(len(fix)):
        n = fix[i]['n']
        for j in range(n):
            power = sum([fix[l]['n'] for l in range(i, len(fix))]) + sum([wait[l]['n'] for l in range(len(wait))]) -j-1
            p121 = p121 + ((1+b)**(f-i-1) * (1+a)**power)**(1/(eta+1)) * (fix[i]['pi'][j])**(eta/(eta+1))
    p122 = 0
    for i in range(len(wait)):
        n = wait[i]['n']
        for j in range(n):
            power = sum(c[i:]) -j-1
            p122 = p122 + ((1+b)**(f-i-1-len(fix)) * (1+a)**power)**(1/(eta+1)) * (wait[i]['pi'][j])**(eta/(eta+1))
    p221 = 0
    for i in range(len(fix)):
        power = sum([fix[l]['n'] for l in range(i, len(fix))]) + sum([wait[l]['n'] for l in range(len(wait))])
        p221 = p221 + ((1+b)**(f-i-1) * (1+a)**power)**(1/(eta+1)) * (fix[i]['omega']*fix[i]['s'])**(eta/(eta+1))
    p222 = 0
    for i in range(len(wait)):
        power = sum(c[i:])
        p222 = p222 + ((1+b)**(f-i-1-len(fix)) * (1+a)**power)**(1/(eta+1)) * (wait[i]['omega']*wait[i]['s'])**(eta/(eta+1))
    p = p11*(p121+p122)+p21*(p221+p222)
    return p


def explode(node):  # [fix, wait, LB]
    for i,e in enumerate(node[1]):
        temp = deepcopy(node)
        temp[0].append(e)  # add an item into fixed
        temp[1].pop(i)  # remove an item from waiting
        temp[2] = LB(temp[0], temp[1])  # Low Bound
        yield temp


def brunch_bound(f):
    theHB, buffer = heuristic(f), list()
    best = [theHB['g'], [], theHB['z']]
    # theHB, buffer, best = {'z':263.2699}, list(), None
    buffer.append(deepcopy([[], g, 0]))  # [fix, wait, LB]
    while(len(buffer)):
        buffer = sorted(buffer, key=lambda e: e[2])  # 先按照LB进行排序,下面再按照展开的节点数进行排序
        buffer = sorted(buffer, key=lambda e: len(e[0]), reverse=True)  # 这样一样可以保证展开节点数多的、LB小的排在前面
        cur = buffer.pop(0)
        if len(cur[0]) is 0:  # node of root
            for node in explode(cur):
                if node[2] > theHB['z']: continue  # cut
                buffer.append(node)
        elif len(cur[0]) is f:  # node of leaf
            if cur[2] > theHB['z']: continue
            if best is None:
                best = cur
            elif cur[2] <= best[2]:
                best = cur
        else:  # node of trunk
            for node in explode(cur):
                if node[2] > theHB['z']: continue
                buffer.append(node)
    return best


def C(seq):
    for i in range(f):
        n = seq[i]['n']
        pxi = list()
        for j in range(n):
            d = {
                'p':    seq[i]['p'][j],
                'xi':   seq[i]['xi'][j]
            }
            pxi.append(d)
        pxi = sorted(pxi, key=lambda x:x['p']*x['xi'])
        for j in range(n):
            power = sum([seq[l]['n'] for l in range(i, len(seq))]) -j-1
            u = (eta*alpha*(1+b)**(f-i-1) * (1+a)**power / beta / pxi[j]['xi'])**(1/(eta+1)) * pxi[j]['p']**(eta/(eta+1))
            try:
                seq[i]['uj'].append(u)
            except KeyError:
                seq[i]['uj'] = list()
                seq[i]['uj'].append(u)
        power = sum([seq[l]['n'] for l in range(i, len(seq))])
        u = (eta*alpha*(1+b)**(f-i-1) * (1+a)**power / gamma*seq[i]['omega'])**(1/(eta+1)) * seq[i]['s']**(eta/(eta+1))
        seq[i]['u'] = u
    p11 = 0
    for i in range(f):
        n = seq[i]['n']
        pxi = list()
        for j in range(n):
            d = {
                'p':    seq[i]['p'][j],
                'xi':   seq[i]['xi'][j]
            }
            pxi.append(d)
        pxi = sorted(pxi, key=lambda x: x['p'] * x['xi'])
        for j in range(n):
            power = sum([seq[l]['n'] for l in range(i, len(seq))]) - j-1
            p11 = p11 + (1+b)**(f-i-1) * (1+a)**power * ((pxi[j]['p'])/(seq[i]['uj'][j]))**eta
    p12 = 0
    for i in range(f):
        power = sum([seq[l]['n'] for l in range(i, len(seq))])
        p12 = p12 + (1+b)**(f-i-1) * (1+a)**power * ((seq[i]['s'])/(seq[i]['u']))**eta
    p = p11 + p12
    return p


def algorithm2(f):
    for i in range(f):
        n = g[i]['n']
        g[i]['pi'] = sorted([g[i]['xi'][j] * g[i]['p'][j] for j in range(n)])
    s, l = [i for i in range(f)], list()
    times = factorial(f)
    c = 0
    for r in range(times):
        try:
            s[c], s[c+1] = s[c+1], s[c]
        except IndexError:
            s[c], s[0] = s[0], s[c]
        c = (c+1) % f
        temp = list()
        for i in range(f):
            idx = s[i]
            temp.append(g[idx])
        l.append(Z(f, temp))
    return min(l)


def initParameters(n, f, instance):
    global g
    i, fname = 0, 'data/%d_%d.txt' % (n, f)
    with open(fname, encoding='utf8') as f:
        for line in f:
            if i is instance:
                g = eval(line.strip())
                return True
            i = i + 1


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: python3 main.py n f instance')
        exit(1)
    n, f, instance = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    for i in range(instance):
        initParameters(n, f, i)

        beg = clock()
        result1 = heuristic(f)
        end = clock()
        print('(%d) Algorithm 3 %f' % (i, (end-beg)))

        beg = clock()
        result2 = brunch_bound(f)
        end = clock()
        print('(%d) B&B %f' % (i, (end-beg)))

        pred = C(result1['g'])
        label = C(result2[0])
        print("(%d) error percentage %f" % (i, ((pred - label) / label)))

        beg = clock()
        pred2 = algorithm2(f)
        end = clock()
        print('(%d) Algorithm 2 %f' % (i, (end-beg)))
