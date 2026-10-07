[README_new_abq.md](https://github.com/user-attachments/files/33147834/README_new_abq.md)
# a-numbers of the Serre-type maximal curves of Abdón–Bezerra–Quoos

Numerical verification code for the paper

> M. Shafighi, *The a-number of the Serre-type maximal curves of Abdón–Bezerra–Quoos*.

## The curves and the result

Let *p* be a prime, *q* = *p*<sup>t</sup>, *Q* = *q*², *m* = (*q*<sup>n</sup> + 1)/(*q* + 1), and let *n* ≥ 3 be odd.
The curve *X*(*q*, *n*) is given by

> *y*<sup>Q</sup> − *y* = *x*<sup>m</sup>

over an algebraic closure of the finite field with *p* elements. Theorem 1 of the paper states that its a-number is

> a(*X*(*q*, *n*)) = (*p* − 1) / (12 *p* (*q* + 1)) · [ (3*q*² + *p* − 2) *q*<sup>n</sup> + *q* ( (*p* − 2) *q*² − (*p* + 1) *q* − 3*p* ) ].

## Files

* `abq_anumber.py` — the three computations that are compared:
  * `formula(q, n)`: the closed formula of Theorem 1 (exact rational arithmetic);
  * `lattice_count(q, n)`: the count of basis monomials *x*<sup>i</sup>*y*<sup>j</sup>*dx* with (*i* mod *p*) + (*j* mod *p*) ≤ *p* − 2, as in Proposition 7;
  * `cartier_rank(q, n)`: brute force. It builds the matrix of the Cartier operator on the basis of holomorphic differentials directly from the definition and computes its rank over F<sub>p</sub>; the a-number is the genus minus this rank. It does not use the monomial formula of Lemma 6.
* `run_checks.py` — runs the comparisons quoted in Remark 10 of the paper.
* `brute_output.txt` — output of `python3 run_checks.py brute`, the run reported in the paper.

## Usage

Python 3.8 or newer, standard library only.

```bash
python3 run_checks.py count      # lattice count vs. Theorem 1 (seconds)
python3 run_checks.py brute      # Cartier rank vs. Theorem 1, all q <= 32 with genus < 10^6 (several minutes)
python3 run_checks.py all
```

A single curve:

```python
from abq_anumber import cartier_rank, lattice_count, formula, params

q, n = 4, 5
g = params(q, n)[3]                                  # the genus
print(g - cartier_rank(q, n), lattice_count(q, n), formula(q, n))   # 409 409 409
```

## License

MIT (see `LICENSE`).
