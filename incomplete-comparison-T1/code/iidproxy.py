import numpy as np, sys, time
Q=3329
P={512:(2,10,4,3),768:(3,10,4,2),1024:(4,11,5,2)}  # k,du,dv,eta1 ; eta2=2
def cbd(n,eta,rng): return (rng.integers(0,2,(n,eta)).sum(1)-rng.integers(0,2,(n,eta)).sum(1))
def icd(AtA,Atb,x0,eta,sw=80):
    x=np.clip(np.round(x0),-eta,eta).astype(float); Ax=AtA@x
    for _ in range(sw):
        ch=0
        for i in range(len(x)):
            a=AtA[i,i]
            if a<=0: continue
            xi=np.clip(np.round(-(Ax[i]-a*x[i]-Atb[i])/a),-eta,eta)
            if xi!=x[i]: Ax+=AtA[:,i]*(xi-x[i]); x[i]=xi; ch+=1
        if ch==0: break
    return x
def recover(pset,Tsize,N,seed=0):
    k,du,dv,eta1=P[pset]; kn=k*256; DIM=2*kn; eta2=2
    duerr=round(Q/2**(du+1)); noise=round(Q/2**(dv+1)); rng=np.random.default_rng(seed)
    x=np.concatenate([cbd(kn,eta1,rng),cbd(kn,eta1,rng)]).astype(float)
    AtA=np.zeros((DIM,DIM)); Atb=np.zeros(DIM); M=N*Tsize; done=0; CH=20000
    while done<M:
        m=min(CH,M-done); done+=m
        Rs=-(cbd(m*kn,eta2,rng).reshape(m,kn)+rng.integers(-duerr,duerr+1,(m,kn)))  # s-part: -(e1+Du)
        Re= cbd(m*kn,eta1,rng).reshape(m,kn)                                        # e-part: r
        R=np.concatenate([Rs,Re],axis=1).astype(float)
        b=R@x+rng.integers(-noise,noise+1,m)
        AtA+=R.T@R; Atb+=R.T@b
    xls=np.linalg.solve(AtA+1e-6*np.eye(DIM),Atb)
    xi=icd(AtA,Atb,xls,eta1); s_ok=int((np.clip(np.round(xi),-eta1,eta1)[:kn]==x[:kn]).all())
    lsq=int((np.clip(np.round(xls),-eta1,eta1)==x).sum())
    return DIM,kn,lsq,int((np.clip(np.round(xi),-eta1,eta1)==x).sum()),s_ok
if __name__=="__main__":
    print("[*] iid short-coefficient proxy (reviewer-validated). full = s-part (kn) exact; e=t-As.",flush=True)
    for p in [512,768,1024]:
        k=P[p][0]; twokn=2*k*256
        print("ML-KEM-%d (2kn=%d):"%(p,twokn),flush=True)
        for mult in [15,30,50,80,120,180]:
            N=mult*twokn; t=time.time()
            DIM,kn,lsq,icdc,sok=recover(p,1,N,0)
            print("  N=%7d (%3dx2kn): LS=%4d/%d ICD=%4d/%d  s-full=%s [%.0fs]"%(N,mult,lsq,DIM,icdc,DIM,bool(sok),time.time()-t),flush=True)
            if sok: break
    print("### |T| sweep at ML-KEM-1024 (find s-full session count per |T|) ###",flush=True)
    twokn=2048
    for T in [1,2,4,8,16,32,51]:
        for mult in [40,70,110,160]:
            N=mult*twokn//T; t=time.time()
            DIM,kn,lsq,icdc,sok=recover(1024,T,N,0)
            if sok:
                print("  |T|=%2d: s-full at N=%6d sessions (%3dx2kn/|T|) [%.0fs]"%(T,N,mult,time.time()-t),flush=True); break
        else:
            print("  |T|=%2d: not full by N=%d"%(T,N),flush=True)
    print("### DONE ###",flush=True)
