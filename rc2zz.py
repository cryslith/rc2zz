#!/usr/bin/env python3

from itertools import *
from random import *

# from the itertools recipes
def grouper(iterable, n, *, fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, fillvalue='x') --> ABC DEF Gxx
    # grouper('ABCDEFG', 3, incomplete='strict') --> ABC DEF ValueError
    # grouper('ABCDEFG', 3, incomplete='ignore') --> ABC DEF
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


## counting conjugacy classes

def conjugacy_representative(k, collapse_zero=False):
    if collapse_zero:
        k = tuple(x for x in k if x != 0)
    d = {}
    n = 1
    result = []
    for x in k:
        if x == 0:
            result.append(x)
            continue
        if x not in d:
            d[x] = n
            n += 1
        result.append(d[x])
    return tuple(result)

def count_conjugacy_classes(m, collapse_zero=False):
    s = set()
    for k in product(range(m+1), repeat=m):
        s.add(conjugacy_representative(k, collapse_zero))
    return len(s)


## inverting encryption

def rc2zz(k):
    S = list(range(256))
    j = 0
    for ctr in range(256):
        j = k[ctr%len(k)]
        S[0], S[j] = S[j], S[0]
    return S


def permutation_to_cycles(S):
    cycles = []
    remaining = set(range(len(S)))
    while remaining:
        x = remaining.pop()
        c = [x]
        y = x
        while True:
            y = S[y]
            if y == x:
                break
            remaining.remove(y)
            c.append(y)
        cycles.append(c)
    return cycles

def cycles_to_permutation(cycles):
    S = [None] * (max(x for c in cycles for x in c) + 1)
    for c in cycles:
        for i in range(len(c)):
            S[c[i]] = c[(i+1)%len(c)]
    return S

def valid_splits(k, C):
    return [j for j in range(1, k+1) if k % j == 0 and gcd(k // j, C) == 1]

def best_partition(values, M):
    best = [(0, [])]
    for i in range(1, M+1):
        b = None
        for v in values:
            if i - v < 0 or best[i - v] is None:
                continue
            if b is None or best[i - v][0] + 1 < b[0]:
                b = (best[i - v][0] + 1, best[i - v][1] + [v])
        best.append(b)
    return best[M][1]

def gcd(a, b):
    if a < 0 or b < 0:
        raise ValueError
    while True:
        if a < b:
            a, b = b, a
        if b == 0:
            return a
        a, b = b, a-b

# a single cycle which is the
# kth root of a product of j cycles each with length C
# j must divide k and gcd(k/j, C) = 1
#
# the resulting cycle has length C_ = j * C
def kth_root_group(k, group, C):
    j = len(group)
    C_ = j * C
    new_cycle = [None] * C_
    for (i, c) in enumerate(group):
        for (l, x) in enumerate(c):
            new_cycle[(l * k + i) % C_] = x
    return new_cycle

# find a kth root with the smallest number of nontrivial cycles
def kth_root(k, cycles):
    output_cycles = [c for c in cycles if len(c) == 1]
    cycles = [c for c in cycles if len(c) > 1]
    cycles_by_lengths = {}
    for c in cycles:
        C = len(c)
        if C not in cycles_by_lengths:
            cycles_by_lengths[C] = []
        cycles_by_lengths[C].append(c)
    for (C, cs) in cycles_by_lengths.items():
        valid = valid_splits(k, C)
        # find the best way to break up the cycles we have into groups
        best = best_partition(valid, len(cs))
        for part in best:
            group = cs[:part]
            cs = cs[part:]
            output_cycles.append(kth_root_group(k, group, C))
    return output_cycles

def find_key(m, S):
    if 256 % m != 0:
        return ValueError('m must divide 256')
    p = 256 // m
    k = find_key_reduced(kth_root(p, permutation_to_cycles(S)))
    if len(k) <= m:
        return k + [0] * (m - len(k))
    raise ValueError(f'no key found (length {len(k)})')

def find_key_reduced(cycles):
    zero_cycle = next(c for c in cycles if 0 in c)
    zpos = zero_cycle.index(0)
    key = zero_cycle[zpos+1:] + zero_cycle[:zpos]
    # skip cycle containing 0 and trivial cycle
    cycles = [c for c in cycles if 0 not in c and len(c) > 1]
    for c in cycles:
        key.extend(c)
    for c in reversed(cycles):
        if len(c) == 1:
            continue
        key.append(c[0])
    key.reverse()
    return key

def rc2zz_reduced(k):
    S = list(range(max(k)+1))
    j = 0
    for ctr in range(len(k)):
        j = k[ctr]
        S[0], S[j] = S[j], S[0]
    return S


def main():
    # for i in range(7):
    #     print((i, count_conjugacy_classes(i)))
    # for i in range(7):
    #     print((i, count_conjugacy_classes(i, True)))

    m = 16
    k = [randrange(256) for _ in range(m)]
    
    print(f'original key: {k}')
    S = rc2zz(k)
    k_found = find_key(m, S)
    S2 = rc2zz(k_found)
    print(f'found key: {k_found}')
    print('key works:', S == S2)


if __name__ == '__main__':
    main()
