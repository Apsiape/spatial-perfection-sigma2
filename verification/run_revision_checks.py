from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,sys,os,json,hashlib,time
V=Path(__file__).resolve().parent
OUT=V/'generated'/'revision-checks';OUT.mkdir(parents=True,exist_ok=True)
shipped={str(p.relative_to(V)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (V/'compiler').glob('*') if p.suffix in ('.zip','.json')}
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
def run(name,args):
 start=time.monotonic();r=subprocess.run([sys.executable,*map(str,args)],cwd=V,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 (OUT/(name+'.log')).write_text(r.stdout,encoding='utf-8',newline='\n')
 result={'name':name,'exit_code':r.returncode,'seconds':round(time.monotonic()-start,3)}
 print(json.dumps(result),flush=True)
 if r.returncode: print(r.stdout[-2000:],flush=True);raise RuntimeError(name)
 return result
jobs=[('good-root',[V/'review/check_good_root.py']),('tree-clock',[V/'review/check_tree_clock.py']),('exact',[V/'compiler/exact_audit.py','--output-dir',OUT/'compiler']),('minimal-compile',[V/'compiler/minimal_compile.py','--output-dir',OUT/'compiler'])]
with ThreadPoolExecutor(max_workers=4) as pool: results=list(pool.map(lambda job:run(*job),jobs))
results.append(run('minimal-structural',[V/'compiler/minimal_structural.py','--output-dir',OUT/'compiler']))
# Exercise the non-Windows branch with WinDLL forbidden, without claiming a native POSIX run.
for script in ['exact_audit','minimal_compile','minimal_structural']:
 portable=OUT/'portable'
 code="import sys,runpy,ctypes; sys.platform='linux'; ctypes.WinDLL=lambda *a,**k: (_ for _ in ()).throw(AssertionError('WinDLL reached')); sys.argv="+repr([str(V/'compiler'/f'{script}.py'),'--output-dir',str(portable)])+"; runpy.run_path(sys.argv[0],run_name='__main__')"
 results.append(run('portable-'+script,['-c',code]))
 target={'exact_audit':'exact_results.json','minimal_compile':'minimal_results.json','minimal_structural':'minimal_structural_results.json'}[script]
 assert json.loads((portable/target).read_text())['peak_working_set_bytes'] is None
 reject=subprocess.run([sys.executable,str(V/'compiler'/f'{script}.py'),'--output-dir',str(V/'compiler')],cwd=V,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 assert reject.returncode==2 and 'shipped artifacts' in reject.stdout
# Verify native and simulated-portable mathematical results against the preserved bundle.
def stable(x):
 if isinstance(x,dict):return {k:stable(v) for k,v in x.items() if k not in ('elapsed_seconds','peak_working_set_bytes')}
 if isinstance(x,list):return [stable(v) for v in x]
 return x
for name in ['exact_results.json','minimal_results.json','minimal_structural_results.json']:
 baseline=json.loads((V/'compiler'/name).read_text())
 for sub in ['compiler','portable']:assert stable(json.loads((OUT/sub/name).read_text()))==stable(baseline),(sub,name)
for rel,digest in shipped.items():assert hashlib.sha256((V/rel).read_bytes()).hexdigest()==digest,rel
summary={'runs':results,'portable_test':'simulated sys.platform=linux, WinDLL forbidden; not a native Linux/macOS execution','bundled_artifacts_unchanged':len(shipped),'mathematical_json_equal_to_shipped':True,'output_dir_rejections':3}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')
print('PASS: all checks, portable branch, artifact immutability, and mathematical comparisons.',flush=True)
