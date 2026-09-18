# Reproduce

## Requirements

- An x86-64 CPU with AVX2 (`grep avx2 /proc/cpuinfo`).
- Rust (stable) and NASM (`nasm --version`); the crate's AVX2 assembly is built via `nasm-rs`.

## 1. Obtain the vulnerable crate

The bug is unmodified upstream code in `pqc_kyber` 0.7.1
([Argyle-Software/kyber](https://github.com/Argyle-Software/kyber)). The only change applied is a
non-default `poc` feature that re-exports a minimal, backend-agnostic surface so an external
harness can act as the attacker (public hashes only) and establish ground truth (`indcpa_dec` with
the real secret key). Enabling the feature does not alter any shipped code path.

```sh
git clone https://github.com/Argyle-Software/kyber rust-kyber
cd rust-kyber
git checkout v0.7.1        # or the 0.7.1 crate release
git apply ../code/rust-kyber-poc.patch
cd ..
```

Place the harness next to it so the layout is:

```
rust-kyber/      (patched crate above)
kyber_poc/
  src/main.rs    (copy of code/main.rs)
  Cargo.toml     (copy of code/Cargo.toml; depends on ../rust-kyber)
```

```sh
mkdir -p kyber_poc/src
cp code/main.rs   kyber_poc/src/main.rs
cp code/Cargo.toml kyber_poc/Cargo.toml
```

## 2. Run

```sh
cd kyber_poc
export RUSTFLAGS="-C target-feature=+aes,+avx2,+sse2,+sse4.1,+bmi2,+popcnt"

# ML-KEM-768 (default), vulnerable AVX2 backend: oracle demo + full key recovery.
cargo run --release --features avx2

# ML-KEM-512 and ML-KEM-1024:
cargo run --release --features avx2,kyber512
cargo run --release --features avx2,kyber1024

# Reference backend (no avx2 feature): implicit rejection works, no oracle, recovery halts.
cargo run --release
```

Set `KYBER_SEED=0x...` to choose the deterministic key (any 64-bit hex). The driver
`code/reproduce.sh` runs 10 distinct keys per parameter set and prints the results table; copy it
next to `kyber_poc/` (or run it from there).

## Expected result (avx2 build)

- **Oracle.** A corrupted ciphertext decapsulates to a different, message-dependent value on the
  avx2 backend than on the reference backend; the attacker recomputes that value from the decrypted
  message using only public data. Plaintext-checking oracle confirmed.
- **Recovery.** The chosen-ciphertext attack recovers the full secret key using only the public
  key and the oracle: 0 coefficient mismatches against ground truth, and the recovered key
  reproduces real decryption on 64 / 64 random ciphertexts. Mean queries per parameter set:

  | Parameter set | Mean queries | Success |
  |---|---|---|
  | ML-KEM-512  | 2876 | 10/10 |
  | ML-KEM-768  | 4263 | 10/10 |
  | ML-KEM-1024 | 9644 | 10/10 |

Captured runs are in [`results/`](results/).
