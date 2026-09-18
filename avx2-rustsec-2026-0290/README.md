# AVX2 `cmov` FO-bypass: key recovery ([RUSTSEC-2026-0290](https://rustsec.org/advisories/RUSTSEC-2026-0290.html) / [-0288](https://rustsec.org/advisories/RUSTSEC-2026-0288.html))

Backend-specific detail for the skipped rejection move in `pqc_kyber` 0.7.1 (and the
`cosmian_kyber` fork). This directory holds the runnable proof of concept, the recovery math, and
the captured transcripts.

## The defect

ML-KEM decapsulation ends with a constant-time conditional move: if re-encryption of the decrypted
message does not match the input ciphertext, the pre-key is overwritten with the secret reject
value `z`, so the returned shared secret `KDF(z || H(c))` is independent of the decrypted message.
The AVX2 backend implements that move as a vectorized blend whose selector mask is a no-op, so the
overwrite never happens. On every ciphertext the pre-key is left equal to `G(m' || h(pk))` with
`m' = Dec(dk, c)`, and decapsulation returns `KDF(G(m' || h(pk))[..32] || H(c))`, a deterministic
function of the decrypted message that the attacker can recompute from public data.

Locus (unmodified upstream code):

- `src/avx2/verify.rs` (the broken `cmov`)
- `src/kem.rs` (decapsulation routing through it)

## Plaintext-checking oracle

The harness demonstrates the oracle before recovering anything. It corrupts a valid ciphertext,
decapsulates it on both backends, and shows:

- the reference backend returns `KDF(z || H(c'))` (implicit rejection works);
- the AVX2 backend returns a different value, and that value equals
  `KDF(G(m' || h(pk))[..32] || H(c'))` recomputed by the attacker from the decrypted message `m'`
  using only the public key;
- a single-bit-wrong message guess produces a different value.

So the AVX2 decapsulation output is a deterministic, attacker-computable function of the decrypted
message: a complete plaintext-checking oracle.

## Recovery

Because the fault accepts every ciphertext, the attacker crafts ciphertexts directly and chooses
the ephemeral vector `u` freely, isolating one secret coefficient per crafted ciphertext and
reading it from the oracle by sweeping the compressed `v` value across a decode boundary. Cost is
`Theta(kn)` queries, with the constant set by the compression precision. The math is in
[`RECOVERY_ANALYSIS.md`](RECOVERY_ANALYSIS.md).

The same harness runs all three parameter sets: the ciphertext bit-packer, the decompressors, the
candidate set `{-eta1, ..., eta1}`, and the probe constant are all derived from the crate
constants, so `d_u = 10 / d_v = 4` (ML-KEM-512, -768) and `d_u = 11 / d_v = 5` (ML-KEM-1024) are
handled by one code path.

## Results

Ten independent keys per parameter set. Each recovered key is checked against the crate's own
ground-truth secret (0 coefficient mismatches) and against reference decryption on 64 random
ciphertexts (64/64).

| Parameter set | Secret dim `kn` | `(d_u, d_v)` | eta1 | Mean queries | Range | Success |
|---|---|---|---|---|---|---|
| ML-KEM-512  | 512  | (10, 4) | 3 | 2876 | 2858-2887 | 10/10 |
| ML-KEM-768  | 768  | (10, 4) | 2 | 4263 | 4213-4318 | 10/10 |
| ML-KEM-1024 | 1024 | (11, 5) | 2 | 9644 | 9550-9784 | 10/10 |

At fixed precision (`d_v = 4`) the query count is linear in the secret dimension (2876 to 4263, a
factor 1.48 against the dimension ratio 1.50). ML-KEM-1024 raises the precision to `d_v = 5`: its
dimension grows only a further factor 1.33, but the query count grows 2.26, the excess coming from
the finer compression grid. Cost is precision-dominated, dimension secondary.

Transcripts:

- [`results/cmov_recovery_sweep.txt`](results/cmov_recovery_sweep.txt): the 10-key sweep above.
- [`results/full_run_mlkem512.txt`](results/full_run_mlkem512.txt),
  [`results/full_run_mlkem768.txt`](results/full_run_mlkem768.txt),
  [`results/full_run_mlkem1024.txt`](results/full_run_mlkem1024.txt): one full verbose run per
  parameter set (oracle proof + recovery + both verifications).

## Reproduce

Build and run steps, including obtaining the vulnerable crate and applying the audit-trail patch,
are in [`REPRODUCE.md`](REPRODUCE.md). In short, on an x86-64 host with AVX2 and NASM:

```sh
cd kyber_poc
export RUSTFLAGS="-C target-feature=+aes,+avx2,+sse2,+sse4.1,+bmi2,+popcnt"
cargo run --release --features avx2               # ML-KEM-768
cargo run --release --features avx2,kyber512      # ML-KEM-512
cargo run --release --features avx2,kyber1024     # ML-KEM-1024
```

## Files

```
README.md               This file.
RECOVERY_ANALYSIS.md    The coefficient-isolation recovery math.
REPRODUCE.md            Build and run instructions.
code/
  main.rs               Deterministic keygen; oracle demonstration; full recovery + verification.
                        Parameter-generic across ML-KEM-512, -768, -1024.
  Cargo.toml            Depends on the patched rust-kyber; --features avx2[,kyber512|,kyber1024].
  rust-kyber-poc.patch  The non-default `poc` feature added to the crate (audit trail; no shipped
                        code path changed). The bug itself is unmodified upstream code.
  reproduce.sh          Driver: 10 keys per parameter set, prints the results table.
results/
  cmov_recovery_sweep.txt   10-key sweep per parameter set.
  full_run_mlkem512.txt     One full verbose run per parameter set.
  full_run_mlkem768.txt
  full_run_mlkem1024.txt
```
