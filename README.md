# a-numbers of the Serre-type maximal curves of Abdón–Bezerra–Quoos

Numerical verification code for the paper

> M. Shafighi, *The a-number of the Serre-type maximal curves of Abdón–Bezerra–Quoos*.

The curves are `X(q,n): y^Q - y = x^m` over `k = \bar F_p`, with `q = p^t`, `Q = q^2`,
`m = (q^n+1)/(q+1)`, `n >= 3` odd. Theorem 1 of the paper gives

    a(X(q,n)) = (p-1) / (12 p (q+1)) * [ (3q^2 + p - 2) q^n + q ((p-2) q^2 - (p+1) q - 3p) ].

## Files

* `abq_anumber.py` — the three computations being compared:
  * `formula(q, n)`: the closed formula of Theorem 1 (exact rational arithmetic);
  * `lattice_count(q, n)`: the count `#{(i,j) in B : (i mod p) + (j mod p) <= p-2}` of Proposition 7;
  * `cartier_rank(q, n)`: brute force. Builds the matrix of the Cartier operator on the basis
    `{x^i y^j dx : (i,j) in B}` of `H^0(X, Omega^1)` directly from the definition and computes
    its rank over `F_p`; `a = g - rank`. It does not use the monomial formula of Lemma 6.
* `run_checks.py` — runs the comparisons quoted in Remark 10.
* `brute_output.txt` — output of `python3 run_checks.py brute` (the run reported in the paper).

## Usage

Python 3.8+, standard library only.

```bash
python3 run_checks.py count      # lattice count vs. Theorem 1 (seconds)
python3 run_checks.py brute      # Cartier rank vs. Theorem 1, all q <= 32 with genus <= 10^6
python3 run_checks.py all
```

A single curve:

```python
from abq_anumber import cartier_rank, lattice_count, formula, params
q, n = 4, 5
g = params(q, n)[3]
print(g - cartier_rank(q, n), lattice_count(q, n), formula(q, n))   # 409 409 409
```

## License

MIT (see `LICENSE`).
