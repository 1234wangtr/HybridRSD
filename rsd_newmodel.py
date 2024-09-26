from math import *
import sys
refresh = True  # refresh:truncation  !refresh:without truncation


def bdt_hybrid(n,k,h):
    return  bdt_hybrid_advanced(n,k,h)

def bdt_hybrid_test():
    '''
        The complexity of different parameters
    '''
    half = False    # half:low weight parameters    !half:large weight parameters
    comp = []

    n, k, h = 2 ** 10, 652, 106
    if half:
        h = 57
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 12, 1589, 172
    if half:
        h = 98
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 14, 3482, 338
    if half:
        h = 198
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 16, 7391, 667
    if half:
        h = 389
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 18, 15336, 1312
    if half:
        h = 760
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 20, 32771, 2467
    if half:
        h = 1419
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 22, 64770, 4788
    if half:
        h = 2735
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    comp.reverse()
    for cp in comp:
        print(cp)


def bdt_hybrid_advanced(n, k, h):
    beta = n // h
    # refresh parameter (with truncation)
    if refresh:
        r = n - k
        n = beta * h
        k = n - r

    print(f"n={n} k={k} h={h} n/h={n / h} beta={beta}")

    best = 1000
    best_u, best_f = -1, -1
    best_u0 = -1
    best_G1,best_G2,best_G3 = 0,0,0
    best_pr = 0
    best_trycomp = 0
    best_n0,best_m0=0,0


    for u in range(beta - 1):
        for f in range(h):
            # \bar{f} each guess \bar{u}+1     h-\bar{f}-g each guess \bar{u}
            for q in range(h-f):    # represents the 'g' in paper
                # f u+1   h-f-q u
                if (n - h - f - (h-q)*u ) <= 0: # guess too many error
                    continue
                log_pr = log2(1 - (u+1) / beta) * f + log2(1 - (u) / beta) * (h-q-f)
                l = (n-k)//(beta-1)

                n0 = k - h - (u+1) * f - (h-f-q)*(u)
                u0 = 0
                G1 = (n-k-l*(beta-u-2))*(n-k)*((n-k-l*(beta-u-2))+n0)
                G1 = max(G1,1)
                log_G1 = log2(G1)
                log_G2 = 0
                log_G3 = 0
                try_comp = 0

                if n0 <= 0: # no need to figure out quadratic equations
                    tot = -log_pr + log_G1
                else:
                    m0 = (h-f-q) * (beta - u - 2) * (beta - u - 1) / 2 + \
                         f * (beta - u - 2) * (beta - u - 3) / 2
                    w = (beta-1) * q
                    u0 = w
                    v0 = n0-w
                    if v0<0:
                        log_G3 = 1000
                        G2 = 2 ** 1000
                        G3 = 2 ** 1000
                    elif v0 <= sqrt(2 * m0) - 2:
                        G2 = v0 ** 4 * n0**2 / 8
                        G2 /= log2(max(v0**2/2,2))
                        G2 = max(1, G2)
                        log_G2 = log2(G2)
                        G3 = 5/4 * (2*v0*v0*w+w*w)+v0**3/log2(max(v0,2))
                        G3 = max(1, G3, v0 * v0 * u0)
                        log_G3 = log2(G3)
                        try_comp = log2(beta) * q
                    else:
                        log_G3 = 1000
                        G2= 2**1000
                        G3 = 2**1000
                    tot = -log_pr + log2(G1+G2+G3*(2**try_comp))
                    #tot = -log_pr + max(log_G1, log_G2, log_G3 + try_comp)
                if tot < best:
                    best = tot
                    best_f = f
                    best_u = u
                    best_u0 = u0
                    best_n0 = n0
                    best_m0 = m0
                    best_G1 = log_G1
                    best_G2 = log_G2
                    best_G3 = log_G3
                    best_trycomp = try_comp
                    best_pr = log_pr

    f = best_f
    u = best_u
    u0 = best_u0
    q = u0//(beta-1)
    n0 = best_n0
    v0 = n0 - u0
    m0 = best_m0
    print(f"!!!best={best} !!!\n f={best_f} u={best_u} u0={best_u0} q={best_u0//(beta-1)}")
    print(f"m0={m0} n0={n0} v0={v0}")
    print(f"pr={best_pr} G1={best_G1} G2={best_G2} G3={best_G3} try_comp={best_trycomp}")
    # print concrete parameter

    print("----------")
    sys.stdout.flush()

    return best


if __name__ == '__main__': 
    bdt_hybrid_test()
    # One can also use bdt_hybrid(n,k,h) to test any parameters
    # n,k,h=10805248, 589760, 1319
    # bdt_hybrid(n,k,h)