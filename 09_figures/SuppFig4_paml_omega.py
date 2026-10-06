"""Fig. S4 - selective constraint on the four island genes.

Bars: omega (dN/dS) under the PAML codeml M0 one-ratio model.
Right-hand column: P value of the branch-model likelihood-ratio test
(two-ratio model with the E1 tip as foreground vs M0, df = 1; "Branch LRT").
It is not the branch-site test.
Run from 09_figures/; writes svgout/Fig. S4.svg.
"""
import os, sys; sys.path.insert(0, '.')
os.makedirs('svgout', exist_ok=True)
from _svglib import *
H=300.0
BLUE='#A6C8E2'; SALMON='#E9A7A3'; LINE='#E08B84'; REDLAB='#B5451B'; TXT='#111111'
# (label, M0 omega, branch-LRT P, highlighted); see 05_selection_pressure/06_run_codeml_four_models.py
# and 07_likelihood_ratio_tests.py
rows=[('GH3_e108',   0.06337,'0.366',False),
      ('MFS',        0.04567,'0.115',False),
      ('GH3_e227',   0.07585,'0.633',True ),
      ('GH1 (bglB)', 0.04247,'0.127',False)]
L,R = 96.0, 560.0          # plot area
T,B = 40.0, 250.0
XMAX=0.15
def xm(v): return L + v/XMAX*(R-L)
XP = xm(0.118)             # P-value column
out=[head(H,'Fig. S4 Selective constraint on the island genes')]
bars=[]; n=len(rows)
step=(B-T)/n
for i,(name,om,p,hl) in enumerate(rows):
    cy=T+step*(i+0.5)
    bh=step*0.46
    bars.append(rect(L, cy-bh/2, xm(om)-L, bh, fill=SALMON if hl else BLUE))
    bars.append(txt(L-8, cy+3.2, name, 9.2, REDLAB if hl else TXT,
                    anchor='end', weight='bold' if hl else 'normal'))
    bars.append(txt(xm(om)+7, cy+3.2, f'ω = {om:.4f}', 9.0, TXT, anchor='start'))
    bars.append(txt(XP, cy+3.2, f'p = {p}', 9.0, TXT, anchor='start', style='italic'))
out.append(g('bars', *bars))
# reference line and column headers
ref=[line(xm(0.1), T-14, xm(0.1), B, LINE, 1.0, dash='4 3'),
     txt(xm(0.1), T-18, 'ω = 0.1', 8.6, LINE),
     txt(XP+4, T-26, 'Branch', 8.6, TXT), txt(XP+4, T-16, 'LRT', 8.6, TXT)]
out.append(g('reference', *ref))
# axis
ax=[line(L, T-4, L, B, TXT, 1.0), line(L, B, R, B, TXT, 1.0)]
for v in (0.00,0.05,0.10,0.15):
    x=xm(v); ax+=[line(x,B,x,B+4,TXT,1.0), txt(x,B+15,f'{v:.2f}',9.0,TXT)]
ax.append(txt((L+R)/2, B+34, 'ω (dN/dS), M0 one-ratio model', 9.6, TXT))
out.append(g('axis', *ax))
out.append(tail())
open('svgout/Fig. S4.svg','w',encoding='utf-8').write(''.join(out))
print('S4 ok')
