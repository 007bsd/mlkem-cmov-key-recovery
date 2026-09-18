#!/usr/bin/env bash
# Reproducible driver for the pqc_kyber AVX2 cmov FO-bypass key recovery.
# Runs 10 distinct deterministic keys per ML-KEM parameter set on the vulnerable
# avx2 backend, checks each recovery against ground truth, and prints the table.
#
# Layout expected (see REPRODUCE.md):
#   rust-kyber/     patched pqc_kyber 0.7.1 (non-default `poc` feature)
#   kyber_poc/      this harness (src/main.rs = code/main.rs, Cargo.toml = code/Cargo.toml)
# Run from the kyber_poc/ directory, or adjust POC_DIR below.
set -euo pipefail
POC_DIR="${POC_DIR:-.}"
cd "$POC_DIR"

# AVX2 intrinsics + NASM path. Host must have AVX2.
export RUSTFLAGS="-C target-feature=+aes,+avx2,+sse2,+sse4.1,+bmi2,+popcnt"

for cfg in "512:avx2,kyber512" "768:avx2" "1024:avx2,kyber1024"; do
  name="${cfg%%:*}"; feat="${cfg##*:}"
  echo "=== ML-KEM-$name  (features: $feat) ==="
  cargo build --release --features "$feat" >/dev/null 2>&1
  tot=0; n=0; ok=0
  for i in $(seq 1 10); do
    seed=$(printf '0x%016X' $(( (0x1000 + i*0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF )))
    out=$(KYBER_SEED="$seed" ./target/release/kyber_poc 2>&1)
    q=$(echo "$out" | grep -oE "decapsulate calls\) : [0-9]+" | grep -oE "[0-9]+$")
    if echo "$out" | grep -q "COMPLETE and VERIFIED"; then res=VERIFIED; ok=$((ok+1)); else res=FAIL; fi
    printf "  key %2d  seed=%s  queries=%s  %s\n" "$i" "$seed" "$q" "$res"
    tot=$((tot+q)); n=$((n+1))
  done
  echo "  mean queries = $((tot/n))   success = $ok/$n"
  echo
done
