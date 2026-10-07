"""
Numerical checks for the a-number of the Serre-type maximal curves
X(q,n): y^Q - y = x^m  (Q = q^2, m = (q^n+1)/(q+1), n >= 3 odd)
of Abdon-Bezerra-Quoos, over k = algebraic closure of F_p, q = p^t.

Three independent ingredients:

  formula(q, n)       closed formula of Theorem 1 (exact rational arithmetic)
  lattice_count(q, n) a = #{(i,j) in B : (i mod p) + (j mod p) <= p-2}
                      (the count of Proposition 7)
  cartier_rank(q, n)  brute force: the matrix of the Cartier operator on the
                      basis {x^i y^j dx : (i,j) in B} of H^0(X, Omega^1) is
                      computed directly from the definition (p-basis expansion
                      of each differential) and its rank is found by Gaussian
                      elimination over F_p. It does NOT use the monomial
                      formula of Lemma 6 (Lucas' theorem / binomial collapse).

a(X) = g - rank(C).
"""
from fractions import Fraction
from math import comb


def char_p(q):
    """Smallest prime dividing q (the characteristic when q is a prime power)."""
    return next(p for p in range(2, q + 1) if q % p == 0)


def params(q, n):
    assert n >= 3 and n % 2 == 1
    p = char_p(q)
    Q = q * q
    m = (q ** n + 1) // (q + 1)
    g = (q - 1) * (q ** n - q) // 2
    return p, Q, m, g


def formula(q, n):
    """Theorem 1."""
    p = char_p(q)
    return Fraction(p - 1, 12 * p * (q + 1)) * (
        (3 * q * q + p - 2) * q ** n + q * ((p - 2) * q * q - (p + 1) * q - 3 * p)
    )


def basis(q, n):
    """B = {(i,j): i >= 0, 0 <= j <= Q-1, Q*i + j*m <= 2g-2}."""
    p, Q, m, g = params(q, n)
    bound = 2 * g - 2
    B = [(i, j) for j in range(Q) if j * m <= bound
         for i in range((bound - j * m) // Q + 1)]
    assert len(B) == g
    return B


def lattice_count(q, n):
    """#{(i,j) in B : i mod p + j mod p <= p-2}, counted by blocks of p in i."""
    p, Q, m, g = params(q, n)
    bound = 2 * g - 2
    total = 0
    for j in range(Q):
        if j * m > bound:
            break
        lim = p - 2 - j % p            # need  i mod p <= lim
        if lim < 0:
            continue
        full, rem = divmod((bound - j * m) // Q + 1, p)   # i = 0..I
        total += full * (lim + 1) + min(rem, lim + 1)
    return total


def cartier_rank(q, n):
    """Rank of the Cartier operator on H^0(X, Omega^1), by brute force.

    Writing y = w^p - x^m with w = y^(Q/p), we have
        x^i y^j = sum_nu (-1)^(j-nu) C(j,nu) w^nu^p ... x^(i+m(j-nu)),
    and the p-basis expansion in x gives, for each nu with
    i + m(j-nu) = p-1 (mod p),  a term  c_nu * w^nu * x^((N_nu-(p-1))/p) dx
    in C(x^i y^j dx).  Powers of w are reduced to the basis 1, y, ..., y^(Q-1)
    using y^Q = y + x^m.  Elements are dicts {(a,b): coeff} meaning x^a y^b.
    """
    p, Q, m, g = params(q, n)
    B = basis(q, n)
    index = {b: k for k, b in enumerate(B)}
    cache = {}

    def ypow(e):                        # y^e in the basis x^a y^b (b < Q)
        if e < Q:
            return {(0, e): 1}
        if e in cache:
            return cache[e]
        res = {}
        for (a, b), c in ypow(e - Q).items():   # y^e = y^(e-Q) * (y + x^m)
            if b + 1 < Q:
                res[(a, b + 1)] = (res.get((a, b + 1), 0) + c) % p
            else:                                # y^Q = y + x^m
                res[(a, 1)] = (res.get((a, 1), 0) + c) % p
                res[(a + m, 0)] = (res.get((a + m, 0), 0) + c) % p
            res[(a + m, b)] = (res.get((a + m, b), 0) + c) % p
        res = {k: v for k, v in res.items() if v}
        cache[e] = res
        return res

    def cartier(i, j):
        out = {}
        for nu in range(j + 1):
            N = i + m * (j - nu)
            if N % p != p - 1:
                continue
            coef = comb(j, nu) * (-1) ** (j - nu) % p
            if coef == 0:
                continue
            s = (N - (p - 1)) // p
            for (a, b), c in ypow(nu * Q // p).items():
                key = (s + a, b)
                out[key] = (out.get(key, 0) + coef * c) % p
        return {k: v for k, v in out.items() if v}

    rows = []
    for (i, j) in B:
        img = cartier(i, j)
        # C maps H^0(Omega^1) into itself: every image monomial must lie in B
        assert all(k in index for k in img), ("image not holomorphic", (i, j))
        rows.append({index[k]: v for k, v in img.items()})

    pivots, rank = {}, 0                # sparse Gaussian elimination over F_p
    for r in rows:
        while r:
            c = min(r)
            if c in pivots:
                pr = pivots[c]
                f = r[c] * pow(pr[c], -1, p) % p
                for k, v in pr.items():
                    nv = (r.get(k, 0) - f * v) % p
                    if nv:
                        r[k] = nv
                    else:
                        r.pop(k, None)
            else:
                pivots[c] = r
                rank += 1
                break
    return rank
