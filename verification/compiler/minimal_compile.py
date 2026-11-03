"""Compile a minimal four-bridge modular-isometry presentation with an adapted
release build_game.py. Exact symbolic projection checks and literal row-zero
contraction. All outputs in the selected output directory; no producer modifications.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,Counter
from itertools import product
import json,zipfile,struct,importlib.util,time,sys,hashlib
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
sys.dont_write_bytecode=True

SRC=BASE/'toolkit'/'build_game.py'
START=time.monotonic()
def deadline():
    if time.monotonic()-START>120:raise TimeoutError('120 second deadline')
params=['e'+str(i) for i in range(17)]+['act'+str(g) for g in range(1,5)]
reflections={'s_'+p for p in params}
def invtoken(t):return t if t in reflections else t[:-1] if t.endswith('~') else t+'~'
def word(w):
    out=[]
    for t in w:
        if out and invtoken(t)==out[-1]:out.pop()
        else:out.append(t)
    return tuple(out)
def wp(a,b):return word(a+b)
def wadj(a):return word(tuple(invtoken(t) for t in a[::-1]))
def poly(*items):
    p=defaultdict(F)
    for w,c in items:p[word(w)]+=F(c)
    return {w:c for w,c in p.items() if c}
def padd(*ps):
    p=defaultdict(F)
    for q,c in ps:
        for w,v in q.items():p[w]+=v*c
    return {w:v for w,v in p.items() if v}
def pmul(p,q):return poly(*[(wp(w,z),a*b) for w,a in p.items() for z,b in q.items()])
def padj(p):return {wadj(w):v for w,v in p.items()}
one={():F(1)}
ep={p:poly(((),F(1,2)),(('s_'+p,),F(1,2))) for p in params}
a=poly((('ap',),F(1,2)),(('am',),F(1,2)))
inputs=[
('support',padd((a,1),(pmul(ep['e13'],pmul(a,ep['e0'])),-1))),
('source',padd((pmul(padj(a),a),1),(ep['e0'],-1))),
('cover',padd(*[(ep['e'+str(i)],1) for i in range(17)],(one,-1)))]
for g in range(1,5):
    prev=0 if g==1 else 1+4*(g-2);st=1+4*(g-1)
    inputs.append(('active'+str(g),padd((ep['act'+str(g)],1),(ep['e'+str(prev)],-1),*[(ep['e'+str(j)],-1) for j in range(st,st+4)])))
words={()}
for p in inputs:
    for w in p[1]:
        for k in range(1,len(w)+1):words.add(w[:k])
for p in params:words.add(('s_'+p,))
for t in ('ap','am','ap~','am~'):words.add((t,))
words=sorted(words,key=lambda w:(len(w),w))
wid={w:'w'+str(i) for i,w in enumerate(words)}
# A formal 3x3 matrix with entries rational group-algebra words.
def fadd(*args):
    out=defaultdict(F)
    for m,c in args:
        for k,v in m.items():out[k]+=v*c
    return {k:v for k,v in out.items() if v}
def scalar(i,j,p):return {(i,j,w):v for w,v in p.items()}
FI={i:scalar(i,i,one) for i in range(3)}
ID=fadd(*[(FI[i],1) for i in range(3)])
def graph(w,i,j,kind='G'):
    c1,c2,c3=(F(1,2),F(1,2),F(1,2)) if kind=='G' else (F(4,5),F(1,5),F(2,5))
    return fadd((FI[i],c1),(FI[j],c2),(scalar(i,j,{w:F(1)}),c3),(scalar(j,i,{wadj(w):F(1)}),c3))
effects={'F'+str(i+1):FI[i] for i in range(3)}
rows=[]
def row(name,terms):
    terms={k:F(v) for k,v in terms.items() if v}
    rows.append({'name':name,'terms':{k:str(v) for k,v in terms.items()}})
def grow(typ,w,i,j):return typ+str(i+1)+str(j+1)+'_'+wid[w]
row('frame',{'F1':1,'F2':1,'F3':1,'I':-1})
for w in words:
    for i,j in ((0,1),(0,2),(1,2)):
        G,H,KG,KH=[grow(t,w,i,j) for t in ('G','H','KG','KH')]
        effects[G]=graph(w,i,j);effects[H]=graph(w,i,j,'H')
        effects[KG]=fadd((FI[i],1),(FI[j],1),(effects[G],-1))
        effects[KH]=fadd((FI[i],1),(FI[j],1),(effects[H],-1))
        row('graphG_'+G,{G:1,KG:1,'F'+str(i+1):-1,'F'+str(j+1):-1})
        row('graphH_'+G,{H:1,KH:1,'F'+str(i+1):-1,'F'+str(j+1):-1})
        row('angle_'+G,{H:5,'F'+str(i+1):-2,'F'+str(j+1):1,G:-4})
triangles=[]
def triangle(u,v,w,name):
    assert wp(u,v)==w
    D='D_'+name;gs=[grow('G',u,0,1),grow('G',v,1,2),grow('G',w,0,2)]
    effects[D]=fadd(*[(effects[g],F(2,3)) for g in gs],(ID,F(-1,3)))
    row(D,{D:3,**{g:-2 for g in gs},'I':1});triangles.append((u,v,w))
triangle((),(),(),'identity')
for w in words:
    if not w:continue
    triangle(w,(),w,'placeL_'+wid[w]);triangle((),w,w,'placeR_'+wid[w])
    if len(w)>1:triangle(w[:-1],(w[-1],),w,'prefix_'+wid[w])
for t in ('ap','am'):triangle((t,),(t+'~',),(),'inverse_'+t)
# Exact parameter binding, using the actual FMS corner projections.
for p in params:
    effects[p]=fadd(*[(scalar(i,i,ep[p]),1) for i in range(3)])
    for i in range(3):
        E,N='EF'+str(i)+'_'+p,'NF'+str(i)+'_'+p
        effects[E]=scalar(i,i,ep[p]);effects[N]=fadd((FI[i],1),(effects[E],-1))
        row(E,{E:1,N:1,'F'+str(i+1):-1})
    row('E_sum_'+p,{**{'EF'+str(i)+'_'+p:1 for i in range(3)},p:-1})
    row('N_sum_'+p,{**{'NF'+str(i)+'_'+p:1 for i in range(3)},p:1,'I':-1})
    for i,j in ((0,1),(0,2),(1,2)):
        names={k:'J'+k+str(i)+str(j)+'_'+p for k in ('11','10','01','00')}
        e=ep[p];ne=padd((one,1),(e,-1))
        for k in names:
            ev=e if k[0]=='1' else ne
            sign=1 if k[1]=='1' else -1
            effects[names[k]]=fadd((scalar(i,i,ev),F(1,2)),(scalar(j,j,ev),F(1,2)),
                (scalar(i,j,ev),F(sign,2)),(scalar(j,i,ev),F(sign,2)))
        row(names['11']+'a',{names['11']:1,names['10']:1,'EF'+str(i)+'_'+p:-1,'EF'+str(j)+'_'+p:-1})
        row(names['11']+'b',{names['00']:1,names['01']:1,'NF'+str(i)+'_'+p:-1,'NF'+str(j)+'_'+p:-1})
        row(names['11']+'c',{names['11']:1,names['01']:1,grow('G',(),i,j):-1})
        row(names['11']+'d',{grow('G',('s_'+p,),i,j):1,names['11']:-1,names['00']:-1})
for name,p in inputs:
    terms=defaultdict(F)
    for w,c in p.items():
        terms[grow('G',w,0,1)]+=2*c
        terms['F1']-=c;terms['F2']-=c
    row('P_'+name,terms)
# Old question lists and a literal forbidden-event stream for four conditional FMS.
variables=sorted(effects);old={s:[{'id':0,'name':'anchor','kind':'anchor','answers':[0]}] for s in ('A','B')}
oldid={}
for name in variables:
    qi=len(old['A']);oldid[name]=qi
    for s in old:old[s].append({'id':qi,'name':name,'kind':'binary','answers':[0,1]})
fms={};triples=lambda par:[tuple(t) for t in product((1,-1),repeat=3) if t[0]*t[1]*t[2]==par]
csel=[(1,1,1),(1,-1,-1)]
for g in range(1,5):
    for side in old:
        for typ,k in ([('R',i) for i in range(3)]+[('T',i) for i in range(2)]) if side=='A' else ([('C',i) for i in range(3)]+[('U',i) for i in range(2)]):
            ans=['off','f']+list(triples(-1 if typ=='C' and k==2 else 1)) if typ in ('R','C') else ['off',1,-1]+[t for t in triples(1) if t!=csel[k]]
            qi=len(old[side]);old[side].append({'id':qi,'name':f'{g}{typ}{k}','kind':'fms','answers':ans,'g':g,'typ':typ,'k':k})
            fms[side,g,typ,k]=qi
for side in old:
    assert all(q['id']==i and len(q['answers'])<=6 for i,q in enumerate(old[side]))
external=set()
def forbidden(qa,ai,question_b,bi):external.add((qa,ai,question_b,bi))
def bindA(q,pred,p):
    for ai,ans in enumerate(old['A'][q]['answers']):
        for bit in (0,1):
            if bool(pred(ans))!=bool(bit):forbidden(q,ai,oldid[p],bit)
def bindB(q,pred,p):
    for bi,ans in enumerate(old['B'][q]['answers']):
        for bit in (0,1):
            if bool(pred(ans))!=bool(bit):forbidden(oldid[p],bit,q,bi)
for g in range(1,5):
    prev=0 if g==1 else 1+4*(g-2);st=1+4*(g-1)
    aq=[old['A'][fms['A',g,t,k]] for t,k in [('R',i) for i in range(3)]+[('T',i) for i in range(2)]]
    bq=[old['B'][fms['B',g,t,k]] for t,k in [('C',i) for i in range(3)]+[('U',i) for i in range(2)]]
    for q in aq:
        bindA(q['id'],lambda a:a!='off','act'+str(g))
        if q['typ']=='R':bindA(q['id'],lambda a:a=='f','e'+str(prev))
    for q in bq:
        bindB(q['id'],lambda a:a!='off','act'+str(g))
        if q['typ']=='C':bindB(q['id'],lambda a:a=='f','e'+str(prev))
    for j,c in enumerate(triples(1)):
        tq=fms['A',g,'T',1 if j==0 else 0]
        bindA(tq,lambda a,c=c:a==c,'e'+str(st+j))
        bindB(fms['B',g,'C',0],lambda a,c=c:a==c,'e'+str(st+j))
    for qa in aq:
        for question_b in bq:
            for ai,a0 in enumerate(qa['answers']):
                for bi,b0 in enumerate(question_b['answers']):
                    reject=(a0=='off')!=(b0=='off')
                    if a0!='off' and b0!='off':
                        ta,tb,i,j=qa['typ'],question_b['typ'],qa['k'],question_b['k']
                        if ta=='R' and tb=='C':
                            reject|=(a0=='f')!=(b0=='f')
                            if isinstance(a0,tuple) and isinstance(b0,tuple):reject|=a0[j]!=b0[i]
                        if ta=='T' and tb=='C' and j==0:
                            reject|=(b0 not in ('f',csel[i])) if a0 in (1,-1) else a0!=b0
                        if ta=='R' and tb=='U' and i==0:
                            reject|=(a0 not in ('f',csel[j])) if b0 in (1,-1) else a0!=b0
                        if ta=='T' and tb=='U' and a0 in (1,-1) and b0 in (1,-1):reject|=a0*b0!=(-1)**(i*j)
                    if reject:forbidden(qa['id'],ai,question_b['id'],bi)
core={'projection_linear_relations':rows,'parameters':params,'unitary_word_labels':[{'id':wid[w],'word':list(w)} for w in words]}
raw=b''.join(struct.pack('<Qq',((16*qa+ai)<<32)|(16*question_b+bi),-1) for qa,ai,question_b,bi in sorted(external))
with zipfile.ZipFile(OUT/'minimal_input.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
    z.writestr('questions.json',json.dumps(old));z.writestr('core.json',json.dumps(core));z.writestr('coefficients.bin',raw)
producer=SRC.read_text()
assert producer.count('assert n%2==0 and maxd==256')==1
assert producer.count("assert len(allq['A'])==len(allq['B'])==169460")==1
producer=producer.replace('assert n%2==0 and maxd==256','assert n%2==0')
producer=producer.replace("assert len(allq['A'])==len(allq['B'])==169460","assert len(allq['A'])==len(allq['B'])")
adapted=OUT/'adapted_build_game.py';adapted.write_text(producer)
sp=importlib.util.spec_from_file_location('adapted_build_game',adapted);mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)
mod.SOURCE=OUT/'minimal_input.zip'
spec,qs,readcore=mod.compile_game()
with zipfile.ZipFile(OUT/'minimal_game.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
    for name,obj in [('game.json',spec),('questions.json',qs),('core.json',readcore)]:z.writestr(name,json.dumps(obj,separators=(',',':')))
with zipfile.ZipFile(OUT/'minimal_game.zip') as z:
    spec=json.loads(z.read('game.json'));qs=json.loads(z.read('questions.json'))
assert len(spec['old_forbidden'])==len(external)
# Form each *serialized* column sum in the formal group algebra.
sums=[]
for plan in spec['rows']:
    coeff=defaultdict(F)
    for name,sgn in plan['slots']:
        coeff['I']+=sgn if name is None else -sgn
        if name is not None:coeff[name]+=2*sgn
    V=fadd(*[(ID if k=='I' else effects[k],v) for k,v in coeff.items() if v])
    source=rows[len(sums)]
    S=fadd(*[(ID if k=='I' else effects[k],F(v)) for k,v in source['terms'].items()])
    assert V==fadd((S,2*plan['den']))
    sums.append(V)
# Exact monomial unitary completions in Q(i); integer complex numbers suffice.
def cmul(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def cconj(x):return x[0],-x[1]
def cadd(x,y):return x[0]+y[0],x[1]+y[1]
def compose(a,b):
    return tuple(a[0][j] for j in b[0]),tuple(cmul(b[1][k],a[1][j]) for k,j in enumerate(b[0]))
def adj(u):
    p=[0]*len(u[0]);v=[None]*len(p)
    for j,i in enumerate(u[0]):p[i]=j;v[i]=cconj(u[1][j])
    return tuple(p),tuple(v)
evals=[]
for m in range(4):
    deadline()
    basis=[(d,x,y) for d in range(m+1) for x in range(2**d) for y in range(2**d)]
    ix={v:i for i,v in enumerate(basis)};t=len(basis);n=17*t
    weights=[1]+[w for w in (2,4,8,16) for _ in range(4)]
    normal=sum((F(4**d,16**d) for d in range(m+1)),F(0))
    eig=[F(weights[c],121*16**d)/normal for c in range(17) for d,x,y in basis]
    unit=(tuple(range(n)),((1,0),)*n);atoms={}
    for p in params:
        if p.startswith('e'):cs={int(p[1:])}
        else:
            g=int(p[3:]);prev=0 if g==1 else 1+4*(g-2);st=1+4*(g-1);cs={prev,*range(st,st+4)}
        atoms['s_'+p]=(tuple(range(n)),tuple((1 if i//t in cs else -1,0) for i in range(n)))
    partial={j:13*t+ix[d+1,2*x,2*y] for j,(d,x,y) in enumerate(basis) if d<m}
    perm=list(range(n));vp=[(1,0)]*n;vm=vp.copy()
    groups=defaultdict(list)
    for j,h in enumerate(eig):groups[h].append(j)
    for cols in groups.values():
        used={partial[j] for j in cols if j in partial}
        holes=[j for j in cols if j not in partial];free=[j for j in cols if j not in used]
        for j in cols:
            if j in partial:perm[j]=partial[j];assert eig[j]==eig[perm[j]]
        for j,r in zip(holes,free):perm[j]=r;vp[j]=(0,1);vm[j]=(0,-1)
    atoms['ap']=(tuple(perm),tuple(vp));atoms['am']=(tuple(perm),tuple(vm))
    for k in ('ap','am'):atoms[k+'~']=adj(atoms[k])
    for u in atoms.values():assert compose(u,adj(u))==unit
    value={():unit}
    for w in words:
        if w:value[w]=compose(value[w[:-1]],atoms[w[-1]])
    def matrix(formal):
        mat={}
        for (fi,fj,w),co in formal.items():
            if w not in value:
                u=unit
                for tok in w:u=compose(u,atoms[tok])
                value[w]=u
            p,v=value[w]
            for j,i in enumerate(p):
                key=(fi*n+i,fj*n+j);z=(co*v[j][0],co*v[j][1])
                mat[key]=cadd(mat.get(key,(F(0),F(0))),z)
        return {k:v for k,v in mat.items() if any(v)}
    failures=[];total=F(0);expected=F(3,4**(m+1)-1)
    for plan,V in zip(spec['rows'],sums):
        if not V:continue
        M=matrix(V)
        if not M:continue
        failures.append(plan['name']);assert plan['name']=='P_source'
        # Direct Hilbert-Schmidt contraction of the signed column sum.
        event=sum((eig[j%n]/3*(v[0]**2+v[1]**2) for (i,j),v in M.items()),F(0))/plan['d']**3
        S2=event*F(plan['d']**3,4*plan['den']**2)
        assert S2==F(2,363)*expected
        assert event==F(8*plan['den']**2,363*plan['d']**3)*expected
        total+=event
        # Verify the actual row-zero game cell is rejecting.
        assert mod.forbidden(spec,qs,plan['question'],0,spec['globals']['X'],0)
    assert failures==['P_source']
    Q=spec['questions_per_party']
    evals.append({'m':m,'tree_dimension':t,'base_dimension':n,'old_dimension':3*n,'compiled_dimension':3*n*spec['scalar_dimension'],
        'failed_rows':failures,'source_L':next(p['den'] for p in spec['rows'] if p['name']=='P_source'),
        'source_d':next(p['d'] for p in spec['rows'] if p['name']=='P_source'),
        'top_mass':str(expected),'sum_rowzero_loss':str(total),'uniform_game_loss':str(total/F(Q*Q)),
        'other_tests':'structural identities proved in report; old FMS event probabilities checked separately'})
    print(json.dumps(evals[-1]),flush=True)
result={'presentation':'four calibrated bridges; one supported modular isometry; no extra transported labels',
 'input_polynomials':len(inputs),'words':len(words),'projection_variables':len(variables),'projection_rows':len(rows),
 'questions_per_party':spec['questions_per_party'],'old_questions':spec['old_questions'],'max_answers':spec['max_answers'],
 'scalar_dimension':spec['scalar_dimension'],'row_dimensions':dict(Counter(p['d'] for p in spec['rows'])),
 'old_forbidden':len(spec['old_forbidden']),'bindings':len(spec['bindings']),'source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),
 'producer_edits':['remove hardcoded maxd==256','remove hardcoded question count==169460'],
 'results':evals,'elapsed_seconds':time.monotonic()-START}

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
(OUT/'minimal_results.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k!='results'},indent=2),flush=True)

