import numpy as np
from kyber_py.ml_kem import ML_KEM_512, ML_KEM_768, ML_KEM_1024
Q=3329; Nn=256; P={512:ML_KEM_512,768:ML_KEM_768,1024:ML_KEM_1024}
def ctr(x): return ((int(x)+Q//2)%Q)-Q//2
def clp(p): return np.array([ctr(c) for c in p.coeffs])
def clv(v,K): return np.array([[ctr(v[i,0].coeffs[a]) for a in range(Nn)] for i in range(K)])
def C1(x):  # Compress_1: bit is 1 iff x mod q in (q/4,3q/4)
    x%=Q; return 1 if (Q/4 < x < 3*Q/4) else 0

def val(pset,nS=30):
    M=P[pset]; K,du,dv=M.k,M.du,M.dv; c1b=K*du*Nn//8; L=1<<dv
    ek,dk=M._k_pke_keygen(M.random_bytes(32))
    s=M.M.decode_vector(dk,K,12,is_ntt=True).from_ntt(); s_c=clv(s,K)
    rho=ek[-32:]; t_hat=M.M.decode_vector(ek[:-32],K,12,is_ntt=True)
    A_T=M._generate_matrix_from_seed(rho,transpose=True)
    match=0; tot=0
    for _ in range(nS):
        m0=M.random_bytes(32); r=M.random_bytes(32); i=Nn-1
        mi=(m0[i//8]>>(i%8))&1
        c=M._k_pke_encrypt(ek,m0,r)
        up=M.M.decode_vector(c[:c1b],K,du).decompress(du); vp=M.R.decode(c[c1b:],dv).decompress(dv)
        upc=clv(up,K)
        sTu=sum(s_c[j][a]*upc[j][(i-a)%Nn]*(1 if a<=i else -1) for j in range(K) for a in range(Nn))
        codes=[int(x)%L for x in M.R.decode(c[c1b:],dv).coeffs]
        for g in range(L):
            deco=round(Q*g/L)%Q
            # harness read: does message bit i survive after overwriting v-coord i with level g?
            pred_accept = (C1(deco - sTu)==mi)
            cc=codes[:]; cc[i]=g
            # real reference decapsulation on the same overwritten ciphertext
            real_accept = (M._k_pke_decrypt(dk,c[:c1b]+M.R(cc).encode(dv))==m0)
            tot+=1; match+= (pred_accept==real_accept)
    return match,tot
if __name__=="__main__":
    import sys; p=int(sys.argv[1]); m,t=val(p)
    print("ML-KEM-%d: harness decode vs real _k_pke_decrypt accept/reject: %d/%d match"%(p,m,t))
