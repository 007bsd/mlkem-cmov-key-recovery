# Full key recovery from a single unverified coordinate (|T| = 1)

This directory demonstrates the paper's central claim: an incomplete Fujisaki-Okamoto ciphertext
comparison that leaves **even one** `v`-coordinate unverified recovers the entire ML-KEM secret
key from a reused key. It is the `|T| = 1` extreme of the incomplete-comparison recovery of
[CVE-2026-10097 / ePrint 2026/1682](https://eprint.iacr.org/2026/1682), which used the full
`|T| ~ 51`-coordinate tail.

## The recovery (decryption-noise regression)

For a valid ciphertext `c = Enc(m0)` the attacker knows `m0`, hence the encryption randomness
`r, e1, e2` and the decompression errors `Du = u' - u`, `Dv = v' - v`. Sweeping the unverified
coordinate `i` over the compression grid and watching the accept/reject flip measures the
decryption noise there to `± q/2^(dv+1)`. Substituting `t = A s + e` into decryption gives a
**known short-coefficient linear equation** in the secret `(s, e)`:

```
delta_i - e2_i - Dv_i  =  ( e^T r  -  s^T (e1 + Du) )_i  =  < (s,e), w_i(m0) >
```

where `w_i(m0)` is built from the negacyclically rotated short vectors `e1+Du` and `r`, so its
entries are small (in `[-eta, eta]` plus a decompression error `<= q/2^(du+1)`), **not** uniform
mod `q`. Each message therefore yields one bounded-error observation of `(s, e)` with short known
coefficients. Because the coefficients are short, this is an ordinary over-determined
least-squares problem, not a lattice problem: stacking `~ 2kn` observations and solving over the
rationals, then rounding to integers, returns the exact `(s, e)`, which the public-key relation
`e = t - A s` certifies. This is the inequality-solving recovery of Hermelink-Pessl-Poeppelmann
(INDOCRYPT 2021) and of 2026/1682; the point here is that one leaked coordinate already supplies
its input.

## Results (against the grid-quantized oracle)

Exact full `(s, e)` recovery at `|T| = 1`, ten random keys per parameter set:

| Parameter set | `2kn` | `dv` | Sessions | Success |
|---|---|---|---|---|
| ML-KEM-512  | 1024 | 4 | 1120 | 10/10 |
| ML-KEM-768  | 1536 | 4 | 1680 | 10/10 |
| ML-KEM-1024 | 2048 | 5 | 2220 | 10/10 |

Recovery is exact (all `2kn` coefficients, zero residual), sharp just above `N = 2kn` sessions.
A `|T|` sweep at ML-KEM-1024 confirms the cost law `sessions ~ 2kn/|T|` across `|T| = 1..51`:

```
|T|    1     2    4    8   16   32   51
N   2088  1044  522  261  130   65   43
```

(the `|T| ~ 51` end is the wolfSSL AVX2 coverage). Full log:
[`results/recovery_sweep.txt`](results/recovery_sweep.txt).

## Reproduce

Python 3 with NumPy and [`kyber-py`](https://github.com/GiacomoPope/kyber-py) (pinned 1.2.0):

```sh
pip install numpy kyber-py
python3 code/recover_T.py --pset 768 --T 1 --max 3000 --cps 1000,2000,3000
```

The oracle is simulated in the `kyber-py` reference implementation: each session takes the exact
decryption noise at the unverified coordinate and rounds it to the compression grid (`± q/2^(dv+1)`),
the accept/reject signal a single-coordinate incomplete comparison exposes. `--T k` sweeps `k`
unverified coordinates; `--pset` selects 512/768/1024. Recovery is checked against the generated
secret (exact match) and is independent of the (separate) cmov coefficient-isolation attack in
[`../avx2-cmov/`](../avx2-cmov/).
