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
        A[r] = (A[r]*pow(int(A[r, c]), q-2, q)) % q
        for rr in range(rows):
            if rr != r and A[rr, c] % q: A[rr] = (A[rr]-A[rr, c]*A[r]) % q
        r += 1
        if r == cols: break
    return r
# Fig 1: rank vs number of sessions, |T|=1, one unchecked coordinate, dim=48 (n=16,k=3)
n, k = 16, 3; dim = k*n; rng = np.random.default_rng(7); forms = []
print("% (sessions, rank) for |T|=1, n=16 k=3 dim=48")
pts = []
for s in range(1, dim+6):
    u = rng.integers(0, q, size=(k, n)); forms.append(form(u, n-1, n, k))
    pts.append((s, rank_mod_q(np.array(forms))))
print(" ".join(f"({s},{r})" for s, r in pts))
