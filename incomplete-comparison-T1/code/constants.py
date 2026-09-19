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
    # I(cell): entropy of grid-cell index of N(0,sigma_sig) quantized to width Delta
    sg=sigma_sig; ext=int(6*sg/Delta)+2; edges=(np.arange(-ext,ext+1))*Delta
    cdf=lambda x:0.5*(1+erf(x/(sg*sqrt(2))))
    pr=np.array([cdf(edges[i+1])-cdf(edges[i]) for i in range(len(edges)-1)]); pr=pr[pr>1e-12]
    I=-(pr*np.log2(pr)).sum()
    print(f"ML-KEM-{p}: sigma_sig={sigma_sig:.1f}  sigma_w in [{sw_min:.2f},{sw_avg:.2f}]  "
          f"c in [{c_min:.2f},{c_avg:.2f}]  Delta={Delta:.0f}  R_suf~{Rlo:.2g}-{Rhi:.2g}  I(cell)={I:.2f}b  "
          f"floor H/I={(2*kn* ( -(0.5*eta1/(2*eta1+1)) ) if False else 0):.0f}")

# Least-squares heuristic R ~ kappa*(Delta/sigma_w)^2 vs measured interface counts
print("\n# LS heuristic: R ~ kappa*(Delta/sigma_w)^2")
meas={512:75000,768:90000,1024:28000}
for p,(k,du,dv,eta1) in P.items():
    kn=k*256; Delta=Q/2**dv; var_r=eta1/2.0; var_du=(Q/2**(du+1))**2/3.0; var_s=1.0+var_du
    sw=(0.5*(var_s+var_r))**0.5; x=(Delta/sw)**2
    print(f"ML-KEM-{p}: (Delta/sw)^2={x:.0f}  measured={meas[p]}  kappa={meas[p]/x:.2f}")
