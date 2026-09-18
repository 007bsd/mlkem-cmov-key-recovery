# Full key recovery from a skipped Fujisaki-Okamoto rejection move in `pqc_kyber` ML-KEM

[RUSTSEC-2026-0290](https://rustsec.org/advisories/RUSTSEC-2026-0290.html) (`pqc_kyber`) and
[RUSTSEC-2026-0288](https://rustsec.org/advisories/RUSTSEC-2026-0288.html) (`cosmian_kyber`, a
stale fork). Both crates are unmaintained and were additionally flagged as such
([RUSTSEC-2026-0289](https://rustsec.org/advisories/RUSTSEC-2026-0289.html),
[RUSTSEC-2026-0287](https://rustsec.org/advisories/RUSTSEC-2026-0287.html)). Advisory PRs:
[Argyle-Software/kyber#121](https://github.com/Argyle-Software/kyber/pull/121),
[Cosmian/kyber#6](https://github.com/Cosmian/kyber/pull/6).

This repository is the reproducibility companion to a forthcoming paper that uses this fault as its
full-coverage instance. The reference will be added here on publication.

## Introduction

ML-KEM (Kyber, [FIPS 203](https://csrc.nist.gov/pubs/fips/203/final)) resists
[chosen-ciphertext attacks](https://en.wikipedia.org/wiki/Chosen-ciphertext_attack) because of one
check inside decapsulation: the receiver re-encrypts the message it just decrypted and returns the
real shared secret only if the result matches the input ciphertext exactly. Any mismatch triggers
implicit rejection, and a pseudorandom value is returned instead. That check is the
[Fujisaki-Okamoto transform](https://eprint.iacr.org/2017/604), the single gate separating the
[IND-CCA2](https://en.wikipedia.org/wiki/Ciphertext_indistinguishability) scheme from the malleable
[IND-CPA](https://en.wikipedia.org/wiki/Ciphertext_indistinguishability) scheme underneath.

The implicit rejection is applied by a constant-time conditional move (`cmov`): on rejection it
overwrites the pre-key with the secret reject value `z`, so the returned shared secret is
independent of the decrypted message. `pqc_kyber`'s hand-written AVX2 backend implemented that
move with a blend whose selector mask is a no-op, so the move never happens. The pre-key is left
equal to `G(m' || h(pk))` with `m' = Dec(dk, c)` on *every* ciphertext, and decapsulation returns
a value that is a deterministic, attacker-computable function of the decrypted message even when
the ciphertext should have been rejected.

That is a complete [plaintext-checking oracle](https://eprint.iacr.org/2019/948). Because the
fault accepts everything, the attacker is not restricted to valid re-encryptions: it crafts
ciphertexts directly, chooses the ephemeral vector `u` freely, isolates one secret coefficient per
crafted ciphertext, and reads it from the oracle. Recovery is `Theta(kn)` decapsulation queries.
The full ML-KEM secret key is recovered end-to-end against the live vulnerable binary at all three
parameter sets (512, 768, 1024), each recovered key verified against the crate's own ground-truth
secret.

The attack needs a reused ML-KEM key: [HPKE](https://www.rfc-editor.org/rfc/rfc9180) recipients,
[KEMTLS](https://eprint.iacr.org/2020/534), or a pinned or embedded key. An ephemeral hybrid
[TLS 1.3](https://www.rfc-editor.org/rfc/rfc8446) key share uses a fresh key per handshake and is
not key-recoverable this way; there the flaw is only a distinguishing break. The advisories are
public, and this is a post-disclosure write-up.

## The fault

One defect, in two crates that share the code, reachable only on the opt-in AVX2 backend.

| | `pqc_kyber` | `cosmian_kyber` |
|---|---|---|
| Advisory | [RUSTSEC-2026-0290](https://rustsec.org/advisories/RUSTSEC-2026-0290.html) | [RUSTSEC-2026-0288](https://rustsec.org/advisories/RUSTSEC-2026-0288.html) |
| Category | crypto-failure (key recovery) | crypto-failure (key recovery) |
| CVSS 3.1 | AV:N/AC:H/PR:N/UI:N/S:U/C:H/I:H/A:N | AV:N/AC:H/PR:N/UI:N/S:U/C:H/I:H/A:N |
| Defect | AVX2 `cmov` blend selector is a no-op; FO rejection copy is skipped | same defect (stale fork of the above) |
| Reachable | `avx2` feature (non-default), x86-64, via public `decapsulate` | same |
| Fixed release | none; crate unmaintained (last 0.7.1) | none; unmaintained fork |
| Also flagged | unmaintained ([RUSTSEC-2026-0289](https://rustsec.org/advisories/RUSTSEC-2026-0289.html)) | unmaintained ([RUSTSEC-2026-0287](https://rustsec.org/advisories/RUSTSEC-2026-0287.html)) |

Full secret-key recovery, 10 independent keys per parameter set, using only the public key and the
decapsulation oracle. Every run recovers the entire key with zero coefficient mismatches against
the ground-truth secret, and the recovered key reproduces reference decryption on 64 of 64 random
ciphertexts.

| Parameter set | Secret dim `kn` | `(d_u, d_v)` | Mean queries | Range | Success |
|---|---|---|---|---|---|
| ML-KEM-512  | 512  | (10, 4) | 2876 | 2858-2887 | 10/10 |
| ML-KEM-768  | 768  | (10, 4) | 4263 | 4213-4318 | 10/10 |
| ML-KEM-1024 | 1024 | (11, 5) | 9644 | 9550-9784 | 10/10 |

Backend-specific detail, math, and reproduction are in
[`avx2-cmov/`](avx2-cmov/).

## FAQ

### What is the vulnerability?

A skipped rejection move in `pqc_kyber`'s AVX2 ML-KEM decapsulation. The constant-time conditional
move that should overwrite the pre-key with the secret reject value on an invalid ciphertext uses a
no-op blend mask, so it never runs. Decapsulation therefore returns a message-dependent value on
every ciphertext, including those a correct implementation rejects, which turns decapsulation into
a complete plaintext-checking oracle and allows full private-key recovery when the key is reused.

### Am I affected?

You are affected if you built `pqc_kyber` (or the `cosmian_kyber` fork) with the non-default
`avx2` feature on an x86-64 target and reuse an ML-KEM private key across decapsulations, as with
HPKE recipients, KEMTLS, or a pinned or embedded key. The default reference backend is not
affected. Because `cargo audit` cannot see feature flags, the advisories flag all affected-version
builds; the exploitable scope is the AVX2 build. Both crates are unmaintained, so there is no fixed
release: migrate to a maintained ML-KEM implementation.

### What can an attacker do?

Recover the entire ML-KEM secret key on any of the three parameter sets, using only the public key
and crafted decapsulation queries, in a few thousand queries (mean 2876 / 4263 / 9644 for
ML-KEM-512 / -768 / -1024). The oracle is a complete plaintext-checking oracle, so the recovery is
exact, not statistical: every one of the ten keys per parameter set was recovered with zero
coefficient mismatches against ground truth.

### How does the attack work?

The `cmov` fault makes decapsulation return a deterministic, public function of the decrypted
message `m' = Dec(dk, c)` on every ciphertext, so the attacker learns each plaintext bit
`b_k = Compress_1(v_k - (s^T u)_k)`. To read the secret coefficient `s_t[i]`, set the ephemeral
vector to a single monomial `u_t = C * X^(n-i)` (all other components zero); the negacyclic product
places `-C * s_t[i]` in message coordinate 0, so `b_0 = Compress_1(v_0 + C * s_t[i])`. Sweeping the
compressed value of `v_0` and watching the oracle bit flip locates a decode boundary, which pins
`C * s_t[i]` to the compression grid and identifies `s_t[i]` in `{-eta1, ..., eta1}`. There are
`kn` coefficients, giving `Theta(kn)` queries, with the constant set by the compression precision.
Full derivation:
[`avx2-cmov/RECOVERY_ANALYSIS.md`](avx2-cmov/RECOVERY_ANALYSIS.md).

### How is this different from the wolfSSL incomplete-comparison key recovery?

The sibling finding in wolfSSL
([CVE-2026-10097](https://www.cve.org/CVERecord?id=CVE-2026-10097), ePrint
[2026/1682](https://eprint.iacr.org/2026/1682),
[007bsd/ml-kem-key-recovery](https://github.com/007bsd/ml-kem-key-recovery)) is a *different* point
in the same class. There the comparison validates all of `u` and skips only a tail of `v`, so the
chosen-`u` ciphertexts a plaintext-checking or key-mismatch attack needs are rejected, and the key
is recovered instead by regressing the decryption noise leaked in the unchecked `v`-coefficients
(partial coverage, roughly 10^5 to 10^6 queries). Here the `cmov` fault validates nothing, `u` is
free, and recovery is the direct coefficient-isolating oracle attack (full coverage, `Theta(kn)`,
a few thousand queries). The two are the full-coverage and partial-coverage extremes of one result:
any rejection-check leak of coverage at least one recovers the whole key. That unification is the
subject of the paper.

### How was it demonstrated?

Against the live pre-fix `pqc_kyber` 0.7.1 AVX2 binary, attacking a reused crate-generated key. The
harness first shows the plaintext-checking oracle: a corrupted ciphertext decapsulates to a
different, message-dependent value on the AVX2 backend than on the reference backend, and the
attacker recomputes that value from the decrypted message using only public data. It then runs the
full recovery and, before reporting any number, checks the recovered key two ways: against the
crate's own ground-truth secret (0 coefficient mismatches) and against reference decryption on 64
random ciphertexts (64/64). This passes on all ten keys at each of ML-KEM-512, -768, and -1024.
Transcripts:
[`avx2-cmov/results/`](avx2-cmov/results/).

### How do I fix it?

Do not use the `avx2` feature of `pqc_kyber` or `cosmian_kyber`. Both crates are unmaintained with
no fixed release, so the durable fix is to migrate to a maintained ML-KEM implementation. The
reference (non-AVX2) backend of these crates is not affected by this defect.

### How serious is it?

High impact, gated preconditions. It is a full private-key recovery (CVSS 3.1 C:H/I:H), but it
needs the opt-in AVX2 build, a reused key, and a plaintext-checking oracle, which is why the vector
carries AC:H. It is not a [side channel](https://en.wikipedia.org/wiki/Side-channel_attack): there
is no timing or power measurement, only a logic bug in a conditional move. It does not affect
correct implementations or the ML-KEM standard itself. Because the fault removes the FO rejection
entirely, the resulting oracle is stronger and the recovery cheaper than the partial-coverage
wolfSSL case above.

## Reproduction

Attack code, math, captured transcripts, and step-by-step reproduction for the key recovery are in
[`avx2-cmov/`](avx2-cmov/). The reduced-parameter model experiments behind the paper's impossibility
result (the spanning lemma, the session-count law, and the rank-climb figure) are in
[`scripts/`](scripts/).

## Credits

The AVX2 `cmov` FO-bypass flaw and the key-recovery demonstration at all three parameter sets are
by [007bsd](https://github.com/007bsd), reported through
[RustSec](https://github.com/rustsec/advisory-db) as
[RUSTSEC-2026-0290](https://rustsec.org/advisories/RUSTSEC-2026-0290.html) (`pqc_kyber`) and
[RUSTSEC-2026-0288](https://rustsec.org/advisories/RUSTSEC-2026-0288.html) (`cosmian_kyber`).

## References

1. RustSec. [RUSTSEC-2026-0290](https://rustsec.org/advisories/RUSTSEC-2026-0290.html) (`pqc_kyber`, key recovery) and [RUSTSEC-2026-0289](https://rustsec.org/advisories/RUSTSEC-2026-0289.html) (unmaintained). Advisory PR [Argyle-Software/kyber#121](https://github.com/Argyle-Software/kyber/pull/121).
2. RustSec. [RUSTSEC-2026-0288](https://rustsec.org/advisories/RUSTSEC-2026-0288.html) (`cosmian_kyber`, key recovery) and [RUSTSEC-2026-0287](https://rustsec.org/advisories/RUSTSEC-2026-0287.html) (unmaintained). Advisory PR [Cosmian/kyber#6](https://github.com/Cosmian/kyber/pull/6).
3. D. Hofheinz, K. Hovelmanns, E. Kiltz. A Modular Analysis of the Fujisaki-Okamoto Transformation. TCC 2017. [ePrint 2017/604](https://eprint.iacr.org/2017/604).
4. P. Ravi, S. S. Roy, A. Chattopadhyay, S. Bhasin. Generic Side-channel Attacks on CCA-secure lattice-based PKE and KEMs. TCHES 2020(3). [ePrint 2019/948](https://eprint.iacr.org/2019/948).
5. Y. Qin, C. Cheng, X. Zhang, Y. Pan, L. Hu, J. Ding. A Systematic Approach and Analysis of Key Mismatch Attacks on Lattice-Based NIST Candidate KEMs. ASIACRYPT 2021. [ePrint 2021/123](https://eprint.iacr.org/2021/123).
6. B. S. Das. Incomplete Ciphertext Comparison in ML-KEM: From an IND-CCA2 Break to Key Recovery. [ePrint 2026/1682](https://eprint.iacr.org/2026/1682). The partial-coverage sibling of this finding.
7. NIST. [FIPS 203: Module-Lattice-Based Key-Encapsulation Mechanism Standard](https://csrc.nist.gov/pubs/fips/203/final).

## Disclosure

Reported through the RustSec advisory database and published as the advisories above. Both crates
are unmaintained with no fixed release; this is a post-disclosure write-up.
