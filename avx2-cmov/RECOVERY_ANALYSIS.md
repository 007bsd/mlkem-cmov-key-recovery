# Coefficient-isolation recovery from the induced oracle

This is the constructive companion to the paper's coefficient-isolation recovery (Appendix A of
"One Unverified Coordinate Is Enough"): it documents exactly what `code/main.rs` computes and why
the query count takes the measured value. Notation follows
ML-KEM / FIPS 203: the ring is `R_q = Z_q[X]/(X^n + 1)` with `n = 256`, `q = 3329`; the secret is
a module vector `s` of rank `k`; `u` is the ephemeral vector and `v` the scalar part of a
ciphertext `c = (u, v)`.

## The linear form read by one query

ML-KEM decryption recovers the message coordinatewise as

```
m' = Compress_1(v - s^T u),
```

where `Compress_1(x) = round(2x/q) mod 2` returns the bit that is 1 exactly when `x mod q` lies in
`(q/4, 3q/4)`. For a crafted ciphertext `c = (u, v)` the k-th recovered bit is

```
b_k = Compress_1( v_k - (s^T u)_k ),      (s^T u)_k = < w_k(u), s >,
```

where `w_k(u)` collects the coefficients of `u` under the negacyclic action `X^n = -1`. Each
query thus reads one adversary-known linear functional of the secret. Because the full-coverage
`cmov` fault accepts every ciphertext, the attacker chooses `u` freely and isolates one secret
coefficient per query (Appendix A), recovering `s` in `Theta(kn)` queries.

Under the full-coverage `cmov` fault the induced oracle returns a known function of
`m' = Dec(dk, c)` on every `c`, so the attacker learns every bit `b_k`. (Under a partial fault it
learns only those `b_k` whose ciphertext lies in the accepted set `S`.)

## Isolating a single secret coefficient

To read the coefficient `s_t[i]` (module component `t`, polynomial index `i`), set the ephemeral
vector to a single monomial:

```
u_t = C * X^(n-i),   u_{t'} = 0 for t' != t,
```

with `C` a decompressed level fixed below. The negacyclic reduction `X^(n-i) * X^i = X^n = -1`
sends the `s_t[i]` term to message coordinate 0:

```
(s^T u)_0 = - C * s_t[i],      so      b_0 = Compress_1( v_0 + C * s_t[i] ).
```

The product also places the other secret coefficients in the remaining coordinates, but
condition 1 below keeps each of those below the decode threshold, so only coordinate 0 responds to
the `v_0` sweep.

Sweeping the decompressed value `v_0` across its `2^{d_v}` compression levels, the oracle bit
`b_0` flips exactly as `v_0 + C * s_t[i]` crosses a decode boundary at distance `q/4`; the flip
position pins `C * s_t[i]` to the grid precision `q / 2^{d_v}`.

Two conditions make the read clean and unique:

1. `eta1 * C < q/4`. This keeps all `O(eta1)` products `C * s_t[i]`, for
   `s_t[i]` in `{-eta1, ..., eta1}`, inside a single decode interval, so the crafted message stays
   in `{0, e_0}` and no other coordinate interferes.
2. `C > q / 2^{d_v}`. This makes consecutive candidate products differ by more than one grid step,
   so the flip identifies `s_t[i]` uniquely.

Both hold at every parameter set. The harness uses `C ~ q/8` when `eta1 = 2` (ML-KEM-768, -1024)
and a smaller `C` when `eta1 = 3` (ML-KEM-512, which has the wider secret range). The probe
u-code `YU` that decompresses to `C` is computed from the crate constants, as is the candidate
set `{-eta1, ..., eta1}`, so a single harness covers all three parameter sets.

## Total cost

Reading one coefficient costs one `v_0`-sweep: at most `2^{d_v}` oracle calls, or `O(d_v)` with a
binary search over the monotone flip. There are `kn` coefficients, giving `Theta(kn)` oracle calls
in total, the constant set by the compression precision `(d_u, d_v)`.

This is the law witnessed by the measurements:

| Parameter set | `kn` | `d_v` | Mean queries | Queries / coeff |
|---|---|---|---|---|
| ML-KEM-512  | 512  | 4 | 2876 | 5.62 |
| ML-KEM-768  | 768  | 4 | 4263 | 5.55 |
| ML-KEM-1024 | 1024 | 5 | 9644 | 9.42 |

At fixed precision (`d_v = 4`) the per-coefficient cost is flat (5.62 vs 5.55), so the total is
linear in `kn`. At ML-KEM-1024 the compression precision rises (`(d_u, d_v)` from `(10, 4)` to
`(11, 5)`) and the per-coefficient cost rises with it (5.55 to 9.42). Cost is precision-dominated,
with the secret dimension a secondary linear factor; we do not claim a closed-form constant.

## What is verified

For each recovered key the harness performs two independent checks:

- **Ground truth.** It compares every recovered coefficient against the real secret obtained from
  the crate's own `indcpa_dec`: 0 mismatches over all `kn` coefficients.
- **Functional.** It re-derives the message for 64 fresh random ciphertexts using only the
  recovered secret and compares against the crate's real decryption: 64 / 64 agree.

Both pass on every one of the 10 keys at each parameter set.
