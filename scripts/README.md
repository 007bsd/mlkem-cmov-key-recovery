# Model experiments

These scripts reproduce the reduced-parameter model experiments behind the paper's impossibility
result. They are independent of the key-recovery harness in [`../avx2-cmov/`](../avx2-cmov/): that
harness recovers a real key from the live vulnerable binary, whereas these witness the underlying
linear-algebra claims at small parameters, where full ML-KEM dimensions would be too large to print
or plot.

The object in all three is the family of linear forms a single leaked rejection coordinate yields.
For a valid ciphertext `c = Enc(m_0)` the ephemeral vector `u = u(m_0)` is a known pseudorandom
vector, and reading coordinate `i` gives the form `(s^T u)_i = <w_i(u), s>`, where `w_i(u)` is a
signed permutation of `u` under the negacyclic action `X^n = -1`. Varying `m_0` rotates `u`, so the
forms collected across sessions probe whether they span the whole secret space.

## Requirements

Python 3 with NumPy. No other dependencies.

```sh
pip install numpy
```

## Scripts

- **`spanning.py`** reproduces the spanning lemma. It collects single-coordinate forms across many
  `m_0` values and computes their rank over `Z_q`, at coverage `|T| = 1` (the worst case, one
  unchecked coordinate) and `|T| = 4`. The forms reach full rank `kn` in every case, so there is no
  provably-safe subspace: any leak of coverage at least one recovers the whole key.

  ```sh
  python3 spanning.py
  ```

- **`sessions.py`** reproduces the session-count law. It measures the average number of sessions
  (distinct `m_0` values) needed to reach full rank, over ten seeds, for several dimensions and
  coverages, and shows it tracks `kn / |T|`.

  ```sh
  python3 sessions.py
  ```

- **`rank_climb.py`** produces the rank-climb figure data: rank against number of sessions at
  coverage `|T| = 1` for `n = 16`, `k = 3` (`dim = 48`). Output is the list of `(sessions, rank)`
  points plotted in the paper.

  ```sh
  python3 rank_climb.py
  ```

Each script prints its results to stdout; the paper's figures and tables use these numbers
directly.
