"""Reproduce the literal m=3,q=3 toy tableau and its LaTeX table."""
from pathlib import Path
import argparse,importlib.util,json,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir','-o',type=Path,default=HERE.parent/'generated/toy')
args=parser.parse_args();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
sp=importlib.util.spec_from_file_location('toy',HERE/'check_tree_clock.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
t=m.tableau(3,3)[3];assert t[3][3]=='f' and t[4][3]=='m'
head=r'''\begin{table}[ht]
\centering\small
\begin{tabular}{@{}cc|cccccccc|c@{}}
\toprule $y$&pulse&0&1&2&3&4&5&6&7&update mask\\\midrule
'''
body=''.join(f'{y}&{int(y%3==0)}&'+'&'.join(row)+f'&{int(y!=7)}'+r'\\'+'\n' for y,row in enumerate(t))
tail=r'''\bottomrule
\end{tabular}
\caption{Depth $m=3$, period $q=3$: the decoder's $8\times8$ toy tableau.
$P_x$ marks column 0, $P_{x,\rm last}$ column 7, $P_y$ row 0 and
$P_{y,\rm last}$ row 7. The last-column neighbour is replaced by blank;
the last-row mask suppresses the toroidal $7\to0$ update.}\label{tab:toy}
\end{table}
'''
(out/'toy-table.tex').write_text(head+body+tail,encoding='utf-8',newline='\n')
(out/'toy-tableau.json').write_text(json.dumps(t,indent=2)+'\n',encoding='utf-8',newline='\n')
print('PASS: pulse at y=3 writes mark at (3,4); generated',out/'toy-table.tex')
