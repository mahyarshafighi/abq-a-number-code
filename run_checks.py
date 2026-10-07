"""
Reproduces the checks quoted in Remark 10 of the paper.

  python3 run_checks.py brute   [--max-genus 1000000]   # slow: minutes
  python3 run_checks.py count                            # fast: seconds
  python3 run_checks.py all
"""
import argparse
import time
from abq_anumber import formula, lattice_count, cartier_rank, params, char_p


def prime_powers(limit):
    out = []
    for q in range(2, limit + 1):
        p = char_p(q)
        r = q
        while r % p == 0:
            r //= p
        if r == 1:
            out.append(q)
    return out


def brute(max_genus, max_q=32):
    bad = total = 0
    for q in prime_powers(max_q):
        for n in range(3, 40, 2):
            g = params(q, n)[3]
            if g > max_genus:
                break
            t = time.time()
            a = g - cartier_rank(q, n)
            f = formula(q, n)
            ok = (a == f)
            total += 1
            bad += not ok
            print(f"q={q:>2} n={n:>2} g={g:>7}  rank-computation a={a:>7}  "
                  f"Theorem 1={f}  {'OK' if ok else 'MISMATCH'}  ({time.time()-t:.1f}s)",
                  flush=True)
    print(f"brute force: {total} pairs, {bad} mismatches")


def count(max_q=127, max_size=10 ** 13):
    bad = total = 0
    for q in prime_powers(max_q):
        for n in range(3, 40, 2):
            if q ** n > max_size:
                break
            total += 1
            if lattice_count(q, n) != formula(q, n):
                bad += 1
                print("MISMATCH", q, n, lattice_count(q, n), formula(q, n))
    print(f"lattice count vs Theorem 1: {total} pairs, {bad} mismatches")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["brute", "count", "all"])
    ap.add_argument("--max-genus", type=int, default=10 ** 6)
    args = ap.parse_args()
    if args.what in ("count", "all"):
        count()
    if args.what in ("brute", "all"):
        brute(args.max_genus)
