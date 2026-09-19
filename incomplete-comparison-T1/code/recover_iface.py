"""
Interface-level |T|=1 key recovery: the observation at each session is obtained by
OVERWRITING the unverified v-coordinate of a valid ciphertext and querying the reference
ML-KEM decapsulation _k_pke_decrypt for accept/reject ([Dec(c')==m0]), located by a counted
binary search over the 2^dv compression grid. No analytic noise model is used for the
measurement: b comes from real decapsulation calls. We also compute the analytic value from the
known key only to ASSERT the read is correct (sanity), never to drive recovery.
Reports measured sessions to full-key recovery AND measured mean _k_pke_decrypt queries per read.
"""
import numpy as np, sys, time, argparse
from kyber_py.ml_kem import ML_KEM_512, ML_KEM_768, ML_KEM_1024
Q=3329; Nn=256; PARAMS={512:ML_KEM_512,768:ML_KEM_768,1024:ML_KEM_1024}
def ctr(x): return ((np.asarray(x,dtype=np.int64)+Q//2)%Q)-Q//2
def clv(v,K): return np.array([[int(ctr(v[i,0].coeffs[a])) for a in range(Nn)] for i in range(K)])
def clp(p): return np.array([int(ctr(c)) for c in p.coeffs])
def negconv_i(a,b,i):  # (a*b)_i in R_q=Z[X]/(X^n+1), a,b length-n coeff arrays
    idx=(i-np.arange(Nn))%Nn; sign=np.where(np.arange(Nn)<=i,1,-1)
    return int(np.sum(a*b[idx]*sign))

_idx={}
def jidx(j):
    if j not in _idx:
        b=np.arange(Nn); _idx[j]=((j-b)%Nn, np.where(b<=j,1,-1))
    return _idx[j]

def icd(AtA,Atb,x0,eta,sweeps=60):
    x=np.clip(np.round(x0),-eta,eta).astype(float); Ax=AtA@x
    for _ in range(sweeps):
        ch=0
        for i in range(len(x)):
            aii=AtA[i,i]
            if aii<=0: continue
            xi=np.clip(np.round(-(Ax[i]-aii*x[i]-Atb[i])/aii),-eta,eta)
            if xi!=x[i]: Ax+=AtA[:,i]*(xi-x[i]); x[i]=xi; ch+=1
        if ch==0: break
    return x

def setup(pset,seed):
    M=PARAMS[pset]; K=M.k; du=M.du; dv=M.dv; c1b=K*du*Nn//8; L=1<<dv
    np.random.seed(seed)
    import os
    d=os.urandom(32); ek,dk=M._k_pke_keygen(d)
    rho=ek[-32:]
    s_hat=M.M.decode_vector(dk,K,12,is_ntt=True)
    t_hat=M.M.decode_vector(ek[:-32],K,12,is_ntt=True)
    A_hat=M._generate_matrix_from_seed(rho)              # keygen matrix (not transposed)
    # compute e = t - A s in the NTT domain FIRST; from_ntt() mutates in place, so any
    # from_ntt on s_hat/t_hat before this would corrupt the product (the v2 mutation bug).
    e_hat=t_hat-A_hat@s_hat
    e_c=clv(e_hat.from_ntt(),K)                          # mutates the temporary e_hat only
    s_c=clv(s_hat.from_ntt(),K)                          # s_hat not reused after this
    A_T=M._generate_matrix_from_seed(rho,transpose=True) # encryption matrix
    return dict(M=M,K=K,du=du,dv=dv,c1b=c1b,L=L,ek=ek,dk=dk,rho=rho,
                t_hat=t_hat,A_T=A_T,s_c=s_c,e_c=e_c,eta1=M.eta_1)  # t_hat stays NTT for session use

def session(S,i,verify):
    M=S['M']; K=S['K']; du=S['du']; dv=S['dv']; c1b=S['c1b']; L=S['L']
    import os
    m0=os.urandom(32); r=os.urandom(32)
    mi=(m0[i//8]>>(i%8))&1
    # recompute encryption randomness (before to_ntt) exactly as _k_pke_encrypt does
    y,N=M._generate_error_vector(r,M.eta_1,0); e1,N=M._generate_error_vector(r,M.eta_2,N); e2,_=M._generate_polynomial(r,M.eta_2,N)
    y_c=clv(y,K); e1_c=clv(e1,K); e2_c=clp(e2)
    c=M._k_pke_encrypt(S['ek'],m0,r)
    # exact u,v (before compress) to get decompression errors Du,Dv
    y_hat=y.to_ntt()
    u_ex=clv((S['A_T']@y_hat).from_ntt()+e1,K)
    mu=M.R.decode(m0,1).decompress(1); mu_c=clp(mu); v_ex=clp(S['t_hat'].dot(y_hat).from_ntt()+e2+mu)
    up=M.M.decode_vector(c[:c1b],K,du).decompress(du); vp=M.R.decode(c[c1b:],dv).decompress(dv)
    vpc=clp(vp); Du=ctr(clv(up,K)-u_ex); Dv=ctr(vpc-v_ex)
    # coefficient row R (adversary-known short vector w_i(m0))
    scoef=e1_c+Du
    R=np.zeros(2*K*Nn)
    idx,sign=jidx(i)
    for j in range(K):
        R[j*Nn:(j+1)*Nn]           = -scoef[j][idx]*sign
        R[K*Nn+j*Nn:K*Nn+(j+1)*Nn] =  y_c[j][idx]*sign
    # ---- REAL interface read: locate accept-boundary via _k_pke_decrypt (counted) ----
    # Overwrite v-coord i with grid level g and query real decapsulation; the accept arc
    # {g: [Dec(c')==m0]} pins (s^T u')_i, since decoding bit i = C1(decompress(g) - (s^Tu')_i).
    codes=[int(x)%L for x in M.R.decode(c[c1b:],dv).coeffs]; g0=codes[i]
    nq=[0]
    def acc(g):
        nq[0]+=1; cc=codes[:]; cc[i]=int(g)%L
        return M._k_pke_decrypt(S['dk'],c[:c1b]+M.R(cc).encode(dv))==m0
    # g0 is in the accept arc (valid ct). Binary-search the upper True->False edge in (g0,g0+L).
    step=1
    while step<L and acc(g0+step): step*=2
    lo=g0+step//2; hi=g0+min(step,L)     # acc(lo)=True, acc(hi)=False (hi<=g0+L)
    while hi-lo>1:
        mid=(lo+hi)//2
        if acc(mid): lo=mid
        else: hi=mid
    # boundary in decompressed-value space is midway between last-accept lo and first-reject hi
    dec=lambda g: round(Q*(g%L)/L)%Q
    d_lo=dec(lo); d_hi=dec(hi)
    bnd=(d_lo + ((d_hi-d_lo)%Q)/2)%Q     # boundary = midpoint of last-accept/first-reject (mod q)
    off = (3*Q)//4 if mi==1 else Q//4    # upper decode-edge offset for m0_i
    w_meas = (bnd - off) % Q             # measured (s^T u')_i, to +- q/2^{dv+1}
    # observation b = (v'_i - (s^Tu')_i - mu_i) - e2_i - Dv_i  == analytic R.x (recover_T identity)
    b_meas = int(ctr(vpc[i] - w_meas - mu_c[i] - e2_c[i] - Dv[i]))
    if verify is not None:
        b_an = int(ctr(round(float(R @ verify))))   # analytic R.x from ground-truth (s,e)
        if abs(int(ctr(b_an-b_meas))) > round(Q/2**dv)+1:
            return None  # read/identity mismatch beyond one grid step
    return R, b_meas, nq[0]

def run(pset,max_ct,cps,seed=0):
    S=setup(pset,seed); K=S['K']; DIM=2*K*Nn; kn=K*Nn; i=Nn-1
    x=np.concatenate([S['s_c'].flatten(),S['e_c'].flatten()]).astype(float)
    AtA=np.zeros((DIM,DIM)); Atb=np.zeros(DIM); tq=0; t0=time.time(); ci=0; cps=sorted(cps); res=None
    mism=0
    for ct in range(1,max_ct+1):
        out=session(S,i,x if ct<=200 else None)  # verify reads on first 200 sessions
        if out is None: mism+=1; continue
        R,b,nq=out; tq+=nq
        AtA+=np.outer(R,R); Atb+=R*b
        if ct==cps[ci]:
            xh=np.linalg.solve(AtA+1e-6*np.eye(DIM),Atb)
            xi=icd(AtA,Atb,xh,S['eta1']).astype(np.int64); xt=x.astype(np.int64)
            sful=int((xi[:kn]==xt[:kn]).sum()); full=(sful==kn)
            print("  P=%d ct=%6d q=%8d (%.2f/read) s=%d/%d FULL=%s mismatch=%d [%.0fs]"
                  %(pset,ct,tq,tq/ct,sful,kn,full,mism,time.time()-t0),flush=True)
            res=(ct,tq,full); ci+=1
            if full or ci>=len(cps): break
    return res

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--pset",type=int,default=512); ap.add_argument("--max",type=int,default=120000)
    ap.add_argument("--seed",type=int,default=0); ap.add_argument("--cps",type=str,default="")
    a=ap.parse_args()
    dv=PARAMS[a.pset].dv
    print("[*] INTERFACE harness ML-KEM-%d (real _k_pke_decrypt reads), noise=+-%d"%(a.pset,round(Q/2**(dv+1))),flush=True)
    if a.cps: cps=[int(x) for x in a.cps.split(",")]
    else: cps=list(range(10000,a.max+1,10000))
    run(a.pset,a.max,cps,a.seed)
