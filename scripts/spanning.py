import numpy as np

# Toy ML-KEM-like ring: Z_q[X]/(X^n+1), module dim k.
# We test the incomplete-comparison fault: the attacker starts from a valid ct=Enc(m0),
# so u = u(m0) is a KNOWN pseudorandom vector; it mutates only UNCHECKED v-coordinates
# (set T) and, per mutated coord i in T, reads the value (s^T u)_i via a threshold
# binary search (a clean linear form in s). Question: do the linear forms collected over
# many m0 (varying u) SPAN the whole secret space, EVEN when |T| is tiny? If yes -> tau=0
# (no provably-safe regime for this fault family).

q = 3329

def negacyclic_form_coeffs(u, i, n, k):
    # coefficient of s_{j,a} in (s^T u)_i  = u_{j,(i-a) mod n} with negacyclic sign
    w = np.zeros(k*n, dtype=np.int64)
    for j in range(k):
        for a in range(n):
            idx = i - a
            if idx >= 0:
                w[j*n + a] = u[j, idx % n] % q
            else:
                w[j*n + a] = (-u[j, (idx + n) % n]) % q
    return w

def rank_mod_q(M, q):
    A = (M.copy() % q).astype(np.int64)
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        piv = None
        for rr in range(r, rows):
            if A[rr, c] % q != 0:
                piv = rr; break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), q-2, q)
        A[r] = (A[r] * inv) % q
        for rr in range(rows):
            if rr != r and A[rr, c] % q != 0:
                A[rr] = (A[rr] - A[rr, c] * A[r]) % q
        r += 1
        if r == cols:
            break
    return r

def experiment(n, k, T, n_m0, seed=0):
    rng = np.random.default_rng(seed)
    dim = k*n
    forms = []
    for _ in range(n_m0):
        u = rng.integers(0, q, size=(k, n))   # pseudorandom u = u(m0)
        for i in T:                            # one readable form per unchecked coord
            forms.append(negacyclic_form_coeffs(u, i, n, k))
    M = np.array(forms, dtype=np.int64)
    return rank_mod_q(M, q), dim, M.shape[0]

for (n,k) in [(16,2),(16,3),(32,2)]:
    dim = k*n
    # WORST case for attacker: only ONE unchecked coordinate (|T|=1)
    T1 = [n-1]
    # need ~dim m0 values to get dim forms
    r1, d1, nf1 = experiment(n,k,T1, n_m0=dim+8, seed=1)
    # coverage |T|=4
    T4 = list(range(n-4,n))
    r4, d4, nf4 = experiment(n,k,T4, n_m0=(dim//4)+4, seed=2)
    print(f"n={n} k={k} dim={dim}")
    print(f"  |T|=1: {nf1} forms from {nf1} m0 -> rank {r1}/{dim}  {'FULL (tau=0)' if r1==dim else 'DEFICIENT (safe subspace!)'}")
    print(f"  |T|=4: {nf4} forms from {nf4//4} m0 -> rank {r4}/{dim}  {'FULL (tau=0)' if r4==dim else 'DEFICIENT'}")
