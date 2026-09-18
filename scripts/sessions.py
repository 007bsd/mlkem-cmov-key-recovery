import numpy as np
q = 3329

def form(u, i, n, k):
    w = np.zeros(k*n, dtype=np.int64)
    for j in range(k):
        for a in range(n):
            idx = i - a
            w[j*n+a] = (u[j, idx % n] if idx >= 0 else -u[j, (idx+n) % n]) % q
    return w

def rank_mod_q(M):
    A = (M % q).astype(np.int64); rows, cols = A.shape; r = 0
    for c in range(cols):
        piv = next((rr for rr in range(r, rows) if A[rr, c] % q), None)
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, c]), q-2, q)) % q
        for rr in range(rows):
            if rr != r and A[rr, c] % q:
                A[rr] = (A[rr] - A[rr, c]*A[r]) % q
        r += 1
        if r == cols: break
    return r

def sessions_to_full(n, k, Tsize, seed):
    rng = np.random.default_rng(seed); dim = k*n
    T = list(range(n-Tsize, n)); forms = []
    for s in range(1, 4*dim):
        u = rng.integers(0, q, size=(k, n))
        for i in T: forms.append(form(u, i, n, k))
        if len(forms) >= dim and rank_mod_q(np.array(forms)) == dim:
            return s
    return -1

print(f"{'n':>3}{'k':>3}{'dim':>5}{'|T|':>5}{'avg_sessions':>14}{'dim/|T|':>9}{'always_full':>13}")
for (n, k) in [(16, 2), (16, 3), (16, 4)]:
    dim = k*n
    for Tsize in [1, 2, 4, 8]:
        if Tsize > n: continue
        res = [sessions_to_full(n, k, Tsize, sd) for sd in range(10)]
        ok = all(r > 0 for r in res)
        avg = np.mean([r for r in res if r > 0]) if ok else float('nan')
        print(f"{n:>3}{k:>3}{dim:>5}{Tsize:>5}{avg:>14.1f}{dim/Tsize:>9.1f}{str(ok):>13}")
