"""
Incomplete-comparison key recovery at arbitrary coverage |T|, all ML-KEM params.
Direct adaptation of the 2026/1682 decryption-noise regression (regress_compressed.py):
verified identity  delta_j - e2_j - c_v_j = (e^T y - s^T(e1+c_u))_j , short known
coefficients, grid measurement to +- q/2^{dv+1}, least-squares recovery of short (s,e).
The only change vs 1682 is that JS (the leaked coordinates) can be a single coordinate,
i.e. |T|=1 -- the case the paper claims but never demonstrated.
"""
import numpy as np, sys, time, argparse
from kyber_py.ml_kem import ML_KEM_512, ML_KEM_768, ML_KEM_1024
Q=3329; Nn=256
PARAMS={512:ML_KEM_512,768:ML_KEM_768,1024:ML_KEM_1024}

def ctr(x): return ((np.asarray(x,dtype=np.int64)+Q//2)%Q)-Q//2

def build(pset):
    M=PARAMS[pset]; K=M.k; DIM=2*K*Nn; du=M.du; dv=M.dv
    c1=K*du*Nn//8; noise=round(Q/2**(dv+1))
    def clv(v): return np.array([[int(ctr(v[i,0].coeffs[a])) for a in range(Nn)] for i in range(K)])
    def clp(p): return np.array([int(ctr(c)) for c in p.coeffs])
    def keygen():
        d=M.random_bytes(32); rho,sigma=M._G(d+bytes([K]))
        s,N=M._generate_error_vector(sigma,M.eta_1,0); e,_=M._generate_error_vector(sigma,M.eta_1,N)
        s_c=clv(s); e_c=clv(e)                       # BEFORE to_ntt (to_ntt mutates in place)
        s1,N1=M._generate_error_vector(sigma,M.eta_1,0); e1k,_=M._generate_error_vector(sigma,M.eta_1,N1)
        t_hat=(M._generate_matrix_from_seed(rho)@s1.to_ntt())+e1k.to_ntt()
        return rho,t_hat,s_c,e_c
    def one_ct(rho,t_hat,_cache={}):
        m=M.random_bytes(32); r=M.random_bytes(32)
        y,N=M._generate_error_vector(r,M.eta_1,0); e1,N=M._generate_error_vector(r,M.eta_2,N); e2,_=M._generate_polynomial(r,M.eta_2,N)
        y_c=clv(y); e1_c=clv(e1); e2_c=clp(e2)       # BEFORE to_ntt
        if _cache.get('rho') is not rho:             # cache the public matrix per key
            _cache['AT']=M._generate_matrix_from_seed(rho,transpose=True); _cache['rho']=rho
        A_T=_cache['AT']; y_hat=y.to_ntt()
        u_exact=(A_T@y_hat).from_ntt()+e1
        mu=M.R.decode(m,1).decompress(1)
        v_exact=t_hat.dot(y_hat).from_ntt()+e2+mu
        u_ex=clv(u_exact); v_ex=clp(v_exact)        # BEFORE compress (compress mutates in place)
        c=u_exact.compress(du).encode(du)+v_exact.compress(dv).encode(dv)
        up=M.M.decode_vector(c[:c1],K,du).decompress(du); vp=M.R.decode(c[c1:],dv).decompress(dv)
        return y_c,e1_c,e2_c, ctr(clv(up)-u_ex), ctr(clp(vp)-v_ex)
    return M,K,DIM,noise,keygen,one_ct

_idx={}
def jidx(j):
    if j not in _idx:
        b=np.arange(Nn); _idx[j]=((j-b)%Nn, np.where(b<=j,1,-1))
    return _idx[j]

def icd(AtA,Atb,x0,eta,sweeps=60):
    # integer coordinate descent over the box [-eta,eta], minimizing ||R x - b||^2
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

def run(pset, Tsize, max_ct, checkpoints, seed=0, verbose=True):
    M,K,DIM,noise,keygen,one_ct=build(pset)
    np.random.seed(seed)
    import kyber_py.ml_kem as _mk  # kyber_py uses os.urandom; reseed via monkeypatch for determinism
    rho,t_hat,s_c,e_c=keygen()
    x=np.concatenate([s_c.flatten(),e_c.flatten()]).astype(np.float64)
    JS=list(range(Nn-Tsize,Nn))                 # |T| unverified coords in the c2 tail
    AtA=np.zeros((DIM,DIM)); Atb=np.zeros(DIM); rng=np.random.default_rng(seed)
    t0=time.time(); ci=0; cps=sorted(checkpoints); res=None
    for ct in range(1,max_ct+1):
        y_c,e1_c,e2_c,c_u,c_v=one_ct(rho,t_hat)
        scoef=e1_c+c_u
        R=np.zeros((len(JS),DIM))
        for ri,j in enumerate(JS):
            idx,sign=jidx(j)
            for i in range(K):
                R[ri, i*Nn:(i+1)*Nn]           = -scoef[i][idx]*sign
                R[ri, K*Nn+i*Nn:K*Nn+(i+1)*Nn] =  y_c[i][idx]*sign
        b=R.dot(x)+rng.integers(-noise,noise+1,size=len(JS))
        AtA+=R.T@R; Atb+=R.T@b
        if ct==cps[ci]:
            xh=np.linalg.solve(AtA+1e-6*np.eye(DIM),Atb)
            eta=M.eta_1; kn=K*Nn; xt=x.astype(np.int64)
            xr=np.clip(np.round(xh),-eta,eta).astype(np.int64); corr=int((xr==xt).sum())
            xi=icd(AtA,Atb,xh,eta).astype(np.int64); corri=int((xi==xt).sum())
            sful=int((xi[:kn]==xt[:kn]).sum())
            full=(sful==kn)                      # s recovered => e=t-As gives the full key
            if verbose:
                print("  P=%d |T|=%d ct=%6d eqns=%7d LS=%d/%d ICD=%d/%d s=%d/%d FULL=%s [%.0fs]"
                      %(pset,Tsize,ct,ct*len(JS),corr,DIM,corri,DIM,sful,kn,full,time.time()-t0)); sys.stdout.flush()
            res=(ct,corri,DIM,full)
            ci+=1
            if full or ci>=len(cps): break
    return res

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--pset",type=int,default=768); ap.add_argument("--T",type=int,default=1)
    ap.add_argument("--max",type=int,default=12000); ap.add_argument("--seed",type=int,default=0)
    ap.add_argument("--cps",type=str,default="2000,4000,6000,8000,10000,12000")
    a=ap.parse_args()
    cps=[int(x) for x in a.cps.split(",")]
    print("[*] incomplete-comparison regression, ML-KEM-%d, |T|=%d, noise=+-%d"%(a.pset,a.T,round(Q/2**(PARAMS[a.pset].dv+1)))); sys.stdout.flush()
    run(a.pset,a.T,a.max,cps,a.seed)
