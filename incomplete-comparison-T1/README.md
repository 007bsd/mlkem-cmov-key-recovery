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

Full key recovery (`s` recovered, all `kn` coefficients; `e = t - A s`), run **end to end against
the reference `_k_pke_decrypt`** (one key per set, each read a binary search over the `2^dv` grid):

| Parameter set | `2kn` | `dv` | Sessions | decaps/read | queries |
|---|---|---|---|---|---|
| ML-KEM-512  | 1024 | 4 | 75,000 | 5.0 | 3.8e5 |
| ML-KEM-768  | 1536 | 4 | 90,000 | 5.0 | 4.5e5 |
| ML-KEM-1024 | 2048 | 5 | 28,000 | 7.0 | 2.0e5 |

The equivalent iid short-coefficient model agrees over five random keys (~70,000 / ~90,000 /
~34,000 sessions) and drives the coverage sweep below.

> Provenance: the interface-harness rows above are `seed 1`, one key per parameter set (a full run
> is ~10^5 sessions, and this machine's memory ceiling made multi-key interface runs impractical, so
> multi-key robustness is shown on the fast model). Only the session count is key-dependent; the
> decaps/read and the bit-for-bit match to `_k_pke_decrypt` are structural.

The right yardstick is not the `2kn` unknown count but the information floor `H(s,e)/I` (secret
entropy over bits fixed per read): about `3.5k / 5.3k / 3.6k` sessions for 512 / 768 / 1024. The
measured counts sit `~20x / 17x / 9.5x` above that floor; the gap is least-squares solver
inefficiency, not missing information (a BP or lattice solver in the HPP21 / Delvaux22 line would
close part of it, not run here). Each read carries well under one bit about the secret (a
±q/2^(dv+1) window on a std-≈47 quantity). ML-KEM-1024 is cheapest: its finer grid (dv=5) halves
the read noise and so raises the per-read information. This is the same regime as 2026/1682's
`|T|~51` recovery (`10^5-10^6` queries). The session count scales as `2kn/|T|`; a coverage sweep at
ML-KEM-1024 confirms it:

```
|T|    1     2    4    8   16   32   51
N   22528 11264 5632 2816 1408 704  441      (one key; sessions to full)
```

Recovery succeeds at every coverage down to one. Full logs:
[`results/recovery_sweep.txt`](results/recovery_sweep.txt).

## The harnesses

- `code/recover_iface.py` is the **interface harness** and produces the table above: every
  observation is the real `_k_pke_decrypt` accept/reject bit `[Dec(c')==m0]` for a ciphertext whose
  unverified `v`-coordinate has been overwritten, located by a counted binary search over the `2^dv`
  grid (measured ~5 decaps/read for dv=4, ~7 for dv=5). It consumes the real decapsulation oracle,
  not a model. Runs in resumable slices (`--state`/`--budget`) since a full run is ~10^5 sessions.
- `code/recover_T.py` runs against the **full kyber-py negacyclic** decapsulation, reading the
  decryption noise at the unverified coordinate to grid precision (LS + integer coordinate descent).
  Its per-coordinate read was validated bit-for-bit against `_k_pke_decrypt` over **every** grid
  level (`480/480`, `480/480`, `960/960` for 512 / 768 / 1024; see `code/ifval2.py`).
- `code/iidproxy.py` runs the **equivalent independent short-coefficient model** (same coefficient
  distribution and grid noise, no per-session crypto). Far faster; it lets us average over five keys
  and sweep coverages, and agrees with the interface harness (70k/90k/34k vs the measured
  75k/90k/28k). The coverage sweep above is from the model.

> Note: an earlier version of this harness read the secret/error coefficients *after* kyber-py's
> in-place `to_ntt()` and `compress()` mutated them, which regressed against a `q`-scale garbage
> functional and made the grid noise negligible (spuriously "exact" recovery at `~1x2kn`). The
> current code reads all short-vector coefficients *before* those calls; the signal std is then
> ≈47 as theory requires. Watch for in-place mutation when adapting.

## Reproduce

Python 3 with NumPy and [`kyber-py`](https://github.com/GiacomoPope/kyber-py) (pinned 1.2.0):

```sh
pip install numpy kyber-py
# interface harness (real _k_pke_decrypt reads); resumable slices, --state saves progress:
python3 code/recover_iface.py --pset 1024 --seed 1 --state st1024.npz --budget 540
python3 code/recover_T.py --pset 1024 --T 1 --max 40000 --cps 30000,40000   # negacyclic (slow)
python3 code/iidproxy.py                                                      # iid model (fast, all params + sweep)
python3 code/ifval2.py 512                                                    # bit-for-bit read validation
```

`--pset` selects 512/768/1024; `recover_T.py --T k` sweeps `k` unverified coordinates. Recovery is
checked against the generated secret; it is independent of the cmov coefficient-isolation attack in
[`../avx2-cmov/`](../avx2-cmov/).
