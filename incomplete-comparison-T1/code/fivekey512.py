import numpy as np
Q=3329; k,du,dv,eta1=2,10,4,3; kn=512; DIM=2*kn; eta2=2
duerr=round(Q/2**(du+1)); noise=round(Q/2**(dv+1))
def cbd(n,eta,rng): return rng.integers(0,2,(n,eta)).sum(1)-rng.integers(0,2,(n,eta)).sum(1)
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
fulls=[]
for seed in range(5):
    rng=np.random.default_rng(100+seed)
    x=np.concatenate([cbd(kn,eta1,rng),cbd(kn,eta1,rng)]).astype(float)
    AtA=np.zeros((DIM,DIM)); Atb=np.zeros(DIM); done=0; full=None
    for N in [55000,65000,75000,85000,95000,110000]:
        while done<N:
            m=min(15000,N-done); done+=m
            Rs=-(cbd(m*kn,eta2,rng).reshape(m,kn)+rng.integers(-duerr,duerr+1,(m,kn)))
            Re=cbd(m*kn,eta1,rng).reshape(m,kn)
            R=np.concatenate([Rs,Re],axis=1).astype(float); b=R@x+rng.integers(-noise,noise+1,m)
            AtA+=R.T@R; Atb+=R.T@b
        xls=np.linalg.solve(AtA+1e-6*np.eye(DIM),Atb); xi=icd(AtA,Atb,xls,eta1)
        if int((np.clip(np.round(xi),-eta1,eta1)!=x).sum())==0: full=N; break
    fulls.append(full if full else 110000)
    print("key %d: full at N=%s"%(seed,fulls[-1]),flush=True)
print("ML-KEM-512 five-key full-recovery N: min=%d max=%d"%(min(fulls),max(fulls)),flush=True)
