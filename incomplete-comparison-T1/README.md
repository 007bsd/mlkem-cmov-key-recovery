# Full key recovery from a single unverified coordinate (|T| = 1)

An incomplete Fujisaki-Okamoto ciphertext comparison in ML-KEM that leaves **even one**
`v`-coordinate unverified recovers the entire secret key from a reused key. This is the `|T| = 1`
extreme of the incomplete-comparison recovery of [CVE-2026-10097 / ePrint
2026/1682](https://eprint.iacr.org/2026/1682), which used the full `|T| ~ 51`-coordinate tail.

## The recovery (decryption-noise regression)

For a valid ciphertext `c = Enc(m0)` the attacker knows `m0`, hence the encryption randomness
`r, e1, e2` and the decompression errors `Du, Dv`. Sweeping the unverified coordinate `i` over
the compression grid and watching the accept/reject flip measures the decryption noise there to
`± q/2^(dv+1)`. Substituting `t = A s + e` gives a **known short-coefficient linear equation** in
the secret `(s, e)`:

```
delta_i - e2_i - Dv_i  =  ( e^T r  -  s^T (e1 + Du) )_i  =  < (s,e), w_i(m0) >
```

The coefficient vector `w_i(m0)` (negacyclically rotated `e1+Du` and `r`) has **short** entries,
so each message yields one bounded-error observation of `(s, e)`. Its error is comparable to its
own size (signal std ≈ 47 against a ± q/2^(dv+1) ≈ ±104 window at dv=4), so this is a noisy
over-determined least-squares problem. We recover the secret `s` by least squares plus integer
coordinate descent in the box `[-eta, eta]`; the error `e` then follows from `e = t - A s`. This
is the inequality-solving recovery of Hermelink-Pessl-Poeppelmann (INDOCRYPT 2021) and of
2026/1682; the point here is that one leaked coordinate already supplies its input.

## Results

Full key recovery (`s` recovered; `e = t - A s`), five random keys per parameter set:

| Parameter set | `2kn` | `dv` | Sessions | ≈ queries | Success |
|---|---|---|---|---|---|
| ML-KEM-512  | 1024 | 4 | ~70,000 | ~2.8e5 | 5/5 |
| ML-KEM-768  | 1536 | 4 | ~90,000 | ~3.6e5 | 5/5 |
| ML-KEM-1024 | 2048 | 5 | ~34,000 | ~1.7e5 | 5/5 |

The counts are a few tens times the `~2kn` information floor, because each read carries well under
one bit about the secret (a ±q/2^(dv+1) window on a std-≈47 quantity). ML-KEM-1024 is cheapest:
its finer grid (dv=5) halves the read noise. This is the same regime as 2026/1682's `|T|~51`
recovery (`10^5-10^6` queries). The session count scales as `2kn/|T|`; a coverage sweep at
ML-KEM-1024 confirms it:

```
|T|    1     2    4    8   16   32   51
N   22528 11264 5632 2816 1408 704  441      (one key; sessions to full)
```

Recovery succeeds at every coverage down to one. Full logs:
[`results/recovery_sweep.txt`](results/recovery_sweep.txt).

## Two harnesses

- `code/recover_T.py` runs the attack against the **full kyber-py negacyclic** decapsulation:
  each session builds a real valid ciphertext, reads the decryption noise at the unverified
  coordinate to grid precision, and accumulates the normal equations (LS + integer coordinate
  descent). Faithful but slow.
- `code/iidproxy.py` runs the **equivalent independent short-coefficient model** (same coefficient
  distribution and grid noise, no per-session crypto). It is far faster, lets us average over keys
  and sweep coverages, and agrees with the negacyclic harness where the two overlap (conditioning,
  not the negacyclic structure, sets the session count). The session-count table above is from the
  model, cross-checked against the negacyclic harness at ML-KEM-512.

> Note: an earlier version of this harness read the secret/error coefficients *after* kyber-py's
> in-place `to_ntt()` and `compress()` mutated them, which regressed against a `q`-scale garbage
> functional and made the grid noise negligible (spuriously "exact" recovery at `~1x2kn`). The
> current code reads all short-vector coefficients *before* those calls; the signal std is then
> ≈47 as theory requires. Watch for in-place mutation when adapting.

## Reproduce

Python 3 with NumPy and [`kyber-py`](https://github.com/GiacomoPope/kyber-py) (pinned 1.2.0):

```sh
pip install numpy kyber-py
python3 code/recover_T.py --pset 1024 --T 1 --max 40000 --cps 30000,40000   # negacyclic (slow)
python3 code/iidproxy.py                                                      # iid model (fast, all params + sweep)
```

`--T k` sweeps `k` unverified coordinates; `--pset` selects 512/768/1024. Recovery is checked
against the generated secret; it is independent of the cmov coefficient-isolation attack in
[`../avx2-cmov/`](../avx2-cmov/).
