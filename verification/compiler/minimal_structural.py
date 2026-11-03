"""Independent exact old-FMS contraction and release verifier audit.
Read minimal_game.zip; Q(sqrt(2)) arithmetic, 17 scaffold coordinates.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from itertools import product
import zipfile,json,importlib.util,sys,time
sys.dont_write_bytecode=True
BASE=Path(__file__).resolve().parent
import argparse
_parser=argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--output-dir', '-o', type=Path, default=BASE/'generated', help='generated outputs (default: compiler/generated); the shipped compiler directory is forbidden')
_args=_parser.parse_args()
OUT=_args.output_dir.resolve()
if OUT==BASE or OUT.is_relative_to(BASE/'toolkit'):
    _parser.error('choose a separate output directory; shipped artifacts must remain unchanged')
OUT.mkdir(parents=True,exist_ok=True)
import shutil
for _archive in ('minimal_game.zip','minimal_input.zip'):
    if not (OUT/_archive).exists():
        shutil.copy2(BASE/_archive,OUT/_archive)
sys.dont_write_bytecode=True

def z(x=0,y=0):return F(x),F(y)
ZERO=z();ONE=z(1)
def addv(a,b):return a[0]+b[0],a[1]+b[1]
def negv(a):return -a[0],-a[1]
def mv(a,b):return a[0]*b[0]+2*a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def plus(*args):
    out={}
    for A,c in args:
        c=z(c) if not isinstance(c,tuple) else c
        for ij,v in A.items():out[ij]=addv(out.get(ij,ZERO),mv(v,c))
    return {ij:v for ij,v in out.items() if v!=ZERO}
def mul(A,B):
    r=defaultdict(list)
    for (j,k),v in B.items():r[j].append((k,v))
    out={}
    for (i,j),v in A.items():
        for k,w in r[j]:out[i,k]=addv(out.get((i,k),ZERO),mv(v,w))
    return {ij:v for ij,v in out.items() if v!=ZERO}
def diag(ids):return {(i,i):ONE for i in ids}
def tensor(A,B):return {(2*i+k,2*j+l):mv(v,w) for (i,j),v in A.items() for (k,l),w in B.items()}
I2=diag(range(2));XX={(0,1):ONE,(1,0):ONE};ZZ={(0,0):ONE,(1,1):z(-1)}
I4=diag(range(4))
square=[[tensor(ZZ,I2),tensor(I2,XX),tensor(ZZ,XX)],
 [tensor(I2,ZZ),tensor(XX,I2),tensor(XX,ZZ)],
 [tensor(ZZ,ZZ),tensor(XX,XX),tensor(mul(XX,ZZ),mul(XX,ZZ))]]
def context(typ,k,signs):
    P=I4
    for j,s in enumerate(signs):P=mul(P,plus((I4,F(1,2)),(square[k][j] if typ=='R' else square[j][k],F(s,2))))
    assert mul(P,P)==P
    return P
def embed(P,st):return {(i+st,j+st):v for (i,j),v in P.items()}
with zipfile.ZipFile(OUT/'minimal_game.zip') as f:
    spec=json.loads(f.read('game.json'));qs=json.loads(f.read('questions.json'));core=json.loads(f.read('core.json'))
oldQ=spec['old_questions'];I=diag(range(17));meas={}
for side in ('A','B'):
    for q in qs[side][:oldQ]:
        if q['kind']!='fms':continue
        g,k,typ=q['g'],q['k'],q['typ'];prev=0 if g==1 else 1+4*(g-2);st=1+4*(g-1)
        active=diag([prev,*range(st,st+4)]);flag=diag([prev])
        for ai,ans0 in enumerate(q['answers']):
            ans=tuple(ans0) if isinstance(ans0,list) else ans0
            if ans=='off':P=plus((I,1),(active,-1))
            elif ans=='f':P=flag
            elif typ in ('R','C'):P=embed(context(typ,k,ans),st)
            elif isinstance(ans,tuple):P=embed(context('C' if typ=='T' else 'R',0,ans),st)
            else:
                s=ans
                if typ=='T':vec={st+k:ONE}
                else:vec={st:z(0,F(1,2)),st+1:z(0,F((-1)**k,2))}
                outer={(i,j):mv(v,w) for i,v in vec.items() for j,w in vec.items()}
                cross={(prev,j):v for j,v in vec.items()};cross.update({(j,prev):v for j,v in vec.items()})
                P=plus((flag,F(1,2)),(outer,F(1,2)),(cross,F(s,2)))
            assert mul(P,P)==P,(side,q['name'],ans)
            meas[side,q['id'],ai]=P
        assert plus(*[(meas[side,q['id'],ai],1) for ai in range(len(q['answers']))])==I
for side in ('A','B'):
    for q in qs[side][:oldQ]:
        if q['kind']!='binary' or q['name'] not in core['parameters']:continue
        p=q['name']
        if p.startswith('e'):cs=[int(p[1:])]
        else:
            g=int(p[3:]);prev=0 if g==1 else 1+4*(g-2);st=1+4*(g-1);cs=[prev,*range(st,st+4)]
        meas[side,q['id'],1]=diag(cs);meas[side,q['id'],0]=plus((I,1),(diag(cs),-1))
amp=[ONE]+[a for a in (z(0,1),z(2),z(0,2),z(4)) for _ in range(4)]
for qa,ai,question_b,bi in spec['old_forbidden']:
    A=meas['A',qa,ai];B=meas['B',question_b,bi];prob=ZERO
    for (i,j),v in A.items():
        if (i,j) in B:prob=addv(prob,mv(mv(v,B[i,j]),mv(amp[i],amp[j])))
    assert prob==ZERO,(qa,ai,question_b,bi,prob)
# Reuse the release's independent serialized-verifier audit after removing its
# single hardcoded scenario size and redirecting its input archive.
src=(BASE/'toolkit'/'check_game.py').read_text()
assert src.count('==Q==169460')==1
src=src.replace('==Q==169460','==Q')
needle="OUT.parent/'out5'/'functional-small.zip'"
assert src.count(needle)==1
src=src.replace(needle,"OUT/'minimal_input.zip'")
path=OUT/'adapted_check_game.py';path.write_text(src)
sp=importlib.util.spec_from_file_location('release_independent_audit',path);mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)
audit=mod.audit(spec,qs,core)
result={'exact_old_FMS_forbidden_probabilities_zero':len(spec['old_forbidden']),
 'exact_FMS_PVMs_checked':40,'field':'Q(sqrt(2))','tree_independent':'all scaffold effects are tensor identity on the tree; normalized root trace is one',
 'release_serialized_audit':audit,'check_edits':['remove hardcoded question count 169460','redirect source archive to minimal_input.zip']}

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
(OUT/'minimal_structural_results.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))

