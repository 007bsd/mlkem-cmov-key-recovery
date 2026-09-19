"""Honest constants for Theorem 1 / the cost discussion, per parameter set.
The observation is b = <(s,e), w>, w = [ s-block: -(e1+Du) ; e-block: r ], all short.
Var(b) = (kn * eta1/2) * (Var(e1+Du) + Var(r)), giving the signal std sigma_sig.
c = E|<d,w>| / ||d||_2 = sigma_w * sqrt(2/pi), with sigma_w the per-entry std of the block the
difference d touches; we report the range over blocks (the union bound's worst case is the
smallest, the typical case the average). Windowed threshold (direct product, no exceptional event):
   R_suf = (Delta / c) * sqrt(2kn) * ln(8*eta1*kn),  Delta = q/2^dv  (units of Z_q; Delta/sigma_w is dimensionless)
I(cell) = entropy (bits) of the width-Delta grid-cell index of the located value ~ N(0, sigma_sig).
Everything is model-dependent; the constant is not tight. Ships as a reproducibility artifact."""
import numpy as np
from math import erf, log2, sqrt, pi, log
Q=3329
P={512:(2,10,4,3),768:(3,10,4,2),1024:(4,11,5,2)}  # k,du,dv,eta1 ; eta2=2
def var_cbd(eta): return eta/2.0
for p,(k,du,dv,eta1) in P.items():
    kn=k*256; Delta=Q/2**dv; eta2=2
    var_r  = var_cbd(eta1)                       # e-block entry ~ r
    var_du = (Q/2**(du+1))**2/3.0                # decompression residual ~ U[-q/2^{du+1}, +...]
    var_s  = var_cbd(eta2)+var_du                # s-block entry ~ -(e1+Du)
    sigma_sig = sqrt(kn*(eta1/2.0)*(var_s+var_r))
    sw_min=sqrt(min(var_s,var_r)); sw_avg=sqrt(0.5*(var_s+var_r))
    c_min=sw_min*sqrt(2/pi); c_avg=sw_avg*sqrt(2/pi)
    Rlo=(Delta/c_avg)*sqrt(2*kn)*log(8*eta1*kn)   # typical
    Rhi=(Delta/c_min)*sqrt(2*kn)*log(8*eta1*kn)   # worst-block
    # I: information per read = E_offset[ H(cell index of N(0,sigma_sig)) ], offset ~ U[0,Delta).
    # (The uniform in-cell offset must be averaged over; fixing offset=0 over-counts, giving ~1.0.)
    sg=sigma_sig; cdf=lambda x:0.5*(1+erf(x/(sg*sqrt(2)))); Hs=[]
    for u in np.linspace(0,Delta,400,endpoint=False):
        ext=int(6*sg/Delta)+3; edges=(np.arange(-ext,ext+1))*Delta+u
        pr=np.array([cdf(edges[i+1])-cdf(edges[i]) for i in range(len(edges)-1)]); pr=pr[pr>1e-12]
        Hs.append(-(pr*np.log2(pr)).sum())
    I=float(np.mean(Hs))
    print(f"ML-KEM-{p}: sigma_sig={sigma_sig:.1f}  sigma_w in [{sw_min:.2f},{sw_avg:.2f}]  "
          f"c in [{c_min:.2f},{c_avg:.2f}]  Delta={Delta:.0f}  R_suf~{Rlo:.2g}-{Rhi:.2g}  I={I:.2f}b")

# Least-squares heuristic R ~ kappa*(Delta/sigma_w)^2 vs measured interface counts
print("\n# LS heuristic: R ~ kappa*(Delta/sigma_w)^2")
meas={512:75000,768:90000,1024:28000}
for p,(k,du,dv,eta1) in P.items():
    kn=k*256; Delta=Q/2**dv; var_r=eta1/2.0; var_du=(Q/2**(du+1))**2/3.0; var_s=1.0+var_du
    sw=(0.5*(var_s+var_r))**0.5; x=(Delta/sw)**2
    print(f"ML-KEM-{p}: (Delta/sw)^2={x:.0f}  measured={meas[p]}  kappa={meas[p]/x:.2f}")

# Identifiability threshold R* = max(small-d, large-d) union-bound contributions (Appendix B):
#   c = sigma_min/sqrt(3)  (Holder + platykurtic 4th moment; valid at every support, unlike Gaussian)
#   small-d:  R ~ (2*Delta/c) ln(8 eta1 kn) sqrt(2kn)
#   large-d:  R ~ 2kn ln(8 eta1 kn) / ln(1/beta0),  beta0 ~ 0.81 (Paley-Zygmund)
print("\n# Identifiability threshold R* (c=sigma_min/sqrt3, Paley-Zygmund beta0=0.81)")
b0=0.81
for p,(k,du,dv,eta1) in P.items():
    kn=k*256; Delta=Q/2**dv
    var_r=eta1/2.0; var_du=(Q/2**(du+1))**2/3.0; var_s=1.0+var_du
    sig_min=np.sqrt(min(var_r,var_s)); c=sig_min/np.sqrt(3); ln=np.log(8*eta1*kn)
    Rs=(2*Delta/c)*ln*np.sqrt(2*kn); Rl=2*kn*ln/np.log(1/b0); R=max(Rs,Rl)
    m={512:75000,768:90000,1024:28000}[p]
    print(f"ML-KEM-{p}: c={c:.2f}  R_small={Rs:.2g}  R_large={Rl:.2g}  R*={R:.2g}  R*/meas={R/m:.1f}x")
