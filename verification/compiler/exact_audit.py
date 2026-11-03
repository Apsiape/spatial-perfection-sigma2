"""Exact independent Fock/clock audit. Standard library only; sparse Q(zeta_3).
Writes only to the selected output directory. 90 second self deadline, below 6 GB.
"""
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import json, time, sys
BASE=Path(__file__).resolve().parent
import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--output-dir', '-o', type=Path, default=BASE/'generated', help='generated outputs (default: compiler/generated); the shipped compiler directory is forbidden')
_args=_parser.parse_args()
OUT=_args.output_dir.resolve()
if OUT==BASE or OUT.is_relative_to(BASE/'toolkit'):
    _parser.error('choose a separate output directory; shipped artifacts must remain unchanged')
OUT.mkdir(parents=True,exist_ok=True)
sys.dont_write_bytecode=True
START=time.monotonic()
def deadline():
    if time.monotonic()-START>90: raise TimeoutError("90-second audit deadline")
class Z:
    __slots__=('a','b')
    def __init__(self,a=0,b=0): self.a,self.b=F(a),F(b)
    @staticmethod
    def cast(x): return x if isinstance(x,Z) else Z(x)
    def __add__(self,x):
        x=Z.cast(x); return Z(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self): return Z(-self.a,-self.b)
    def __sub__(self,x): return self+-Z.cast(x)
    def __rsub__(self,x): return Z.cast(x)+-self
    def __mul__(self,x):
        x=Z.cast(x)
        return Z(self.a*x.a-self.b*x.b,self.a*x.b+self.b*x.a-self.b*x.b)
    __rmul__=__mul__
    def conj(self): return Z(self.a-self.b,-self.b)
    def norm2(self): return self.a*self.a-self.a*self.b+self.b*self.b
    def __truediv__(self,x):
        x=Z.cast(x); n=x.norm2(); y=self*x.conj()
        return Z(y.a/n,y.b/n)
    def __pow__(self,n):
        if n<0:return (Z(1)/self)**(-n)
        x,y=Z(1),self
        while n:
            if n&1:x=x*y
            y=y*y;n//=2
        return x
    def __bool__(self):return bool(self.a or self.b)
    def __eq__(self,x):
        x=Z.cast(x);return self.a==x.a and self.b==x.b
    def __repr__(self):return str(self.a) if not self.b else f"({self.a})+({self.b})*zeta3"
def clean(a): return {ij:z for ij,z in a.items() if z}
def plus(*args):
    a={}
    for b,c in args:
        for ij,z in b.items(): a[ij]=a.get(ij,Z())+z*c
    return clean(a)
def mul(a,b):
    right=defaultdict(list)
    for (j,k),v in b.items():right[j].append((k,v))
    out={}
    for (i,j),v in a.items():
        for k,w in right[j]:out[i,k]=out.get((i,k),Z())+v*w
    return clean(out)
def star(a):return {(j,i):v.conj() for (i,j),v in a.items()}
def diag(values):return {(i,i):Z.cast(v) for i,v in enumerate(values) if v}
def model(m,q):
    basis=[(d,x,y,j) for d in range(m+1) for x in range(2**d) for y in range(2**d) for j in range(q)]
    ix={w:i for i,w in enumerate(basis)};n=len(basis);I=diag([1]*n)
    zeta=Z(-1) if q==2 else Z(0,1)
    A={}
    for a in (0,1):
        for b in (0,1):
            A[a,b]={(ix[d+1,2*x+a,2*y+b,j],i):Z(1) for i,(d,x,y,j) in enumerate(basis) if d<m}
    R=diag([int(d==0) for d,x,y,j in basis]);P=diag([int(d==m) for d,x,y,j in basis])
    U={(ix[d,(x+1)%(2**d),y,j],i):Z(1) for i,(d,x,y,j) in enumerate(basis)}
    V={(ix[d,x,(y+1)%(2**d),j],i):Z(1) for i,(d,x,y,j) in enumerate(basis)}
    X=diag([int(x==0) for d,x,y,j in basis]);Y=diag([int(y==0) for d,x,y,j in basis])
    W=diag([zeta**y for d,x,y,j in basis]);C=diag([F(1,2**y) for d,x,y,j in basis])
    E=diag([int(y%q==0) for d,x,y,j in basis])
    dv=[Z() if y%q==0 else Z(F(1,2**(y+1)))/(zeta**y-1) for d,x,y,j in basis]
    D=diag(dv);Zop=diag([zeta]*n)
    J=diag([zeta**j for d,x,y,j in basis])
    K={(ix[d,x,y,(j+1)%q],i):Z(1) for i,(d,x,y,j) in enumerate(basis)}
    normal=sum((F(4**d,16**d) for d in range(m+1)),F(0))
    weights=[F(1,q*16**d)/normal for d,x,y,j in basis];h=diag(weights)
    tests=[]
    def eq(name,left,right={}):
        deadline();assert left==right,(m,q,name,len(plus((left,1),(right,-1))))
        tests.append(name)
    for name,T in [('U',U),('V',V),('W',W),('J',J),('K',K),('Z',Zop)]:
        eq(name+' unitary',mul(star(T),T),I);eq(name+' unitary reverse',mul(T,star(T)),I)
    for name,T in [('R',R),('X',X),('Y',Y),('E',E),('top',P)]:
        eq(name+' projection',mul(T,T),T);eq(name+' selfadjoint',star(T),T)
    for name,T in [('U',U),('V',V),('X',X),('Y',Y),('W',W),('C',C),('E',E),('D',D),('J',J),('K',K),('Z',Zop),('R',R)]:
        eq(name+' centralizer',mul(h,T),mul(T,h))
    for alpha,T in A.items():
        eq('top isometry '+str(alpha),mul(star(T),T),plus((I,1),(P,-1)))
        eq('modular '+str(alpha),mul(h,T),plus((mul(T,h),F(1,16))))
        a,b=alpha
        for name,L in [('J',J),('K',K),('Z',Zop)]:
            eq(name+' creation commutator '+str(alpha),mul(L,T),mul(T,L))
        eq('U recursion '+str(alpha),mul(U,T),A[1,b] if a==0 else mul(A[0,b],U))
        eq('V recursion '+str(alpha),mul(V,T),A[a,1] if b==0 else mul(A[a,0],V))
        eq('X recursion '+str(alpha),mul(X,T),mul(T,X) if a==0 else {})
        eq('Y recursion '+str(alpha),mul(Y,T),mul(T,Y) if b==0 else {})
        rhs=mul(T,mul((Zop if b else I),mul(W,W)))
        eq('W recursion '+str(alpha),mul(W,T),rhs)
        eq('C recursion '+str(alpha),mul(C,T),plus((mul(T,mul(C,C)),F(1,2**b))))
    for alpha,T in A.items():
        for beta,S in A.items():
            if alpha!=beta:eq('cross isometry '+str(alpha)+str(beta),mul(star(T),S))
    eq('range sum',plus(*[(mul(T,star(T)),1) for T in A.values()]),plus((I,1),(R,-1)))
    for name,T in [('U',U),('V',V),('X',X),('Y',Y),('W',W),('C',C)]:
        eq(name+' root',mul(T,R),R)
    eq('Weyl',mul(J,K),mul(Zop,mul(K,J)))
    eq('Z J',mul(Zop,J),mul(J,Zop));eq('Z K',mul(Zop,K),mul(K,Zop))
    eq('UV',mul(U,V),mul(V,U))
    eq('reader support',mul(plus((W,1),(I,-1)),E))
    eq('reader inverse',mul(D,plus((W,1),(I,-1))),plus((mul(C,plus((I,1),(E,-1))),F(1,2))))
    eq('last X',mul(star(U),mul(X,U)),diag([int(x==2**d-1) for d,x,y,j in basis]))
    eq('last Y',mul(star(V),mul(Y,V)),diag([int(y==2**d-1) for d,x,y,j in basis]))
    assert all(0<v.a<=1 and not v.b for (i,j),v in C.items()) # injective finite truncation
    assert max((z.norm2() for z in dv),default=F(0))<=1
    eq('reader kernel exactly',E,diag([int(zeta**y==1) for d,x,y,j in basis]))
    mass=sum((weights[i] for i,(d,x,y,j) in enumerate(basis) if d==m),F(0))
    assert mass==F(3,4**(m+1)-1)
    # Explicit clock-sensitive deterministic CA, no halting symbol reached.
    # At epsilon=0 the unique moving front advances; at epsilon=1 it freezes.
    # An origin start s seeds the front at t=1 regardless of epsilon at t=0.
    gamma=['s','b','f','z','H']
    def rule(l,c,r,e):
        if c=='H':return 'H'
        if c=='s':return 'b'
        if c=='z':return 'z'
        if c=='f':return 'z' if e else 'b'
        if l=='s':return 'f'
        if l=='f' and not e:return 'f'
        return 'b'
    tile={}
    for d in range(m+1):
        width=2**d;row=['s']+['b']*(width-1)
        for y in range(width):
            for x,c in enumerate(row):tile[d,x,y]=c
            if y+1<width:
                row=[rule('#' if x==0 else row[x-1],row[x],row[x+1] if x+1<width else 'b',int(y%q==0)) for x in range(width)]
    pc={c:diag([int(tile[d,x,y]==c) for d,x,y,j in basis]) for c in gamma}
    eq('CA PVM',plus(*[(pc[c],1) for c in gamma]),I)
    for c in gamma:eq('CA projection '+c,mul(pc[c],pc[c]),pc[c])
    eq('CA start',mul(pc['s'],mul(X,Y)),mul(X,Y))
    eq('CA blank',mul(pc['b'],mul(plus((I,1),(X,-1)),Y)),mul(plus((I,1),(X,-1)),Y))
    eq('CA halt',pc['H'])
    xlast=mul(star(U),mul(X,U));ylast=mul(star(V),mul(Y,V))
    ls={'#':X};ns={}
    for c in gamma:
        ls[c]=mul(mul(U,mul(pc[c],star(U))),plus((I,1),(X,-1)))
        ns[c]=plus((mul(mul(star(U),mul(pc[c],U)),plus((I,1),(xlast,-1))),1),(xlast,int(c=='b')))
    eps=[plus((I,1),(E,-1)),E];no_last=plus((I,1),(ylast,-1))
    accum={c:{} for c in gamma}
    for l in ['#']+gamma:
        for c in gamma:
            for r in gamma:
                base=mul(ls[l],mul(pc[c],ns[r]))
                if not base:continue
                for e in (0,1):
                    dest=rule(l,c,r,e)
                    accum[dest]=plus((accum[dest],1),(mul(base,eps[e]),1))
    for c in gamma:eq('CA update '+c,mul(mul(star(V),mul(pc[c],V)),no_last),mul(accum[c],no_last))
    # Four explicit shared bridges, copies of noncommuting labels, and exact
    # lifted relation table. The field is unchanged; only sparse dimensions grow.
    chosen=[0,1,5,9,13];scweights=[1]+[w for w in (2,4,8,16) for _ in range(4)]
    def at(T,r,s):return {(r*n+i,s*n+j):v for (i,j),v in T.items()}
    hs=diag([scweights[c]*v/F(121) for c in range(17) for v in weights])
    bridges=[at(I,chosen[g],chosen[g-1]) for g in range(1,5)]
    v=bridges[0]
    for wg in bridges[1:]:v=mul(wg,v)
    eq('chain source',mul(star(v),v),at(I,0,0))
    eq('chain range',mul(v,star(v)),at(I,13,13))
    eq('chain ratio16',mul(hs,v),plus((mul(v,hs),16)))
    for Lname,L in [('U',U),('V',V),('X',X),('Y',Y),('W',W),('C',C),('E',E),('D',D),('J',J),('K',K),('Z',Zop),('R',R)]:
        for g,wg in enumerate(bridges,1):
            eq('shared bridge '+Lname+str(g),mul(at(L,chosen[g],chosen[g]),wg),mul(wg,at(L,chosen[g-1],chosen[g-1])))
        eq('shared chain '+Lname,mul(at(L,13,13),v),mul(v,at(L,0,0)))
    lifted={alpha:at(T,13,0) for alpha,T in A.items()}
    for alpha,T in lifted.items():
        eq('lift centralizer '+str(alpha),mul(hs,T),mul(T,hs))
        eq('extraction '+str(alpha),mul(star(v),T),at(A[alpha],0,0))
        eq('lift source '+str(alpha),mul(star(T),T),at(plus((I,1),(P,-1)),0,0))
    eq('lift range',plus(*[(mul(T,star(T)),1) for T in lifted.values()]),at(plus((I,1),(R,-1)),13,13))
    for b in (0,1):
        eq('translated U zero-bit '+str(b),mul(at(U,13,13),lifted[0,b]),lifted[1,b])
        eq('translated U one-bit '+str(b),mul(at(U,13,13),lifted[1,b]),mul(lifted[0,b],at(U,0,0)))
    # J,K share the bridge even though they do not commute.
    assert mul(J,K)!=mul(K,J)

    return {'m':m,'q':q,'root_dimension':n,'checks':len(tests),'failed_original_rows':['A_%d%d* A_%d%d = I'%(a,b,a,b) for a,b in A],
        'top_mass_per_row':str(mass),'sum_four_squared_source_defects':str(4*mass),
        'D_norm_squared':str(max((z.norm2() for z in dv),default=F(0))),'CA':'explicit nonhalting pulse-sensitive front, not the unspecified universal CA',
        'q_exceeds_all_times':q>2**m-1}
results=[]
for m in range(5):
    for q in (2,3):
        results.append(model(m,q));print(json.dumps(results[-1]),flush=True)
result={'arithmetic':'exact Fraction pairs in Q(zeta_3), with zeta_2=-1','results':results,'bad_roots_at_q_2_or_3':[],
        'elapsed_seconds':time.monotonic()-START}

def peak_bytes():
    if sys.platform != 'win32':
        return None  # Windows working-set telemetry is unavailable.
    import ctypes
    from ctypes import wintypes
    class PMC(ctypes.Structure):
        _fields_=[('cb',wintypes.DWORD),('PageFaultCount',wintypes.DWORD)]+[(x,ctypes.c_size_t) for x in ('PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage')]
    pm=PMC();pm.cb=ctypes.sizeof(pm)
    kernel=ctypes.WinDLL('kernel32',use_last_error=True);psapi=ctypes.WinDLL('psapi',use_last_error=True)
    kernel.GetCurrentProcess.restype=wintypes.HANDLE
    psapi.GetProcessMemoryInfo.argtypes=[wintypes.HANDLE,ctypes.c_void_p,wintypes.DWORD]
    assert psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(),ctypes.byref(pm),pm.cb)
    assert pm.PeakWorkingSetSize<6*1024**3
    return pm.PeakWorkingSetSize

result['peak_working_set_bytes']=peak_bytes()
(OUT/'exact_results.json').write_text(json.dumps(result,indent=2))
print('PASS',len(results),'models',sum(r['checks'] for r in results),'exact checks')


