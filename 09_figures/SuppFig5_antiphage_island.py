"""Fig. S5 - strain-unique anti-phage defence island of S. vesiculosa E1.

Arrow map of the type I R-M locus (hsdR, hsdS, hsdM), the two E1-unique
DUF262/DUF1524 genes, relB, a phage integrase and a WYL-domain regulator, with
Pfam domains below each arrow. Coordinates are the E1 Prokka annotation.
Run from 09_figures/ after `mkdir -p svgout`; writes svgout/Fig. S5.svg.
"""
import sys; sys.path.insert(0,'.')
from _svglib import *
H=182.0
RM='#2C5F8D'; UNIQ='#E9000D'; MOB='#B5B5B5'; MOB2='#CFCFCF'; TA='#74A9CF'; WYL='#9C7BAA'
TXT='#222222'; SUB='#555555'; GREY='#8A8A8A'
RY=72.0; GH=26.0
X0, X1 = 30.0, 470.0
K0, K1 = 2814.8, 2830.2
def km(kb): return X0 + (kb-K0)/(K1-K0)*(X1-X0)
# (name, domain annotation, start kb, end kb, colour, points left, italic name)
genes=[('',              'Phage_integrase',  2815.2, 2816.6, MOB,  True,  False),
       ('',              'DUF45',            2816.9, 2817.5, MOB2, True,  False),
       ('hsdR',          'HSDR_N,ResIII',    2817.6, 2820.4, RM,   True,  True ),
       ('',              'DUF1524,DUF262',   2820.4, 2822.0, UNIQ, True,  False),
       ('',              'DUF262',           2822.0, 2823.5, UNIQ, True,  False),
       ('hsdS',          'Methylase_S',      2823.5, 2825.0, RM,   True,  True ),
       ('hsdM',          'HsdM_N,N6_Mtase',  2825.0, 2827.1, RM,   True,  True ),
       ('relB',          'RelB',             2827.2, 2827.8, TA,   True,  True ),
       ('',              '',                 2827.9, 2828.5, MOB2, True,  False),
       ('',              'WYL',              2828.6, 2829.9, WYL,  False, False)]
out=[head(H,'Fig. S5 Strain-unique anti-phage defence island in S. vesiculosa E1')]
body=[line(X0-4, RY, X1+4, RY, '#E4E4E4', 2.4, cap='round')]
for name, dom, a, b, col, lft, ital in genes:
    xa, xb = km(a), km(b)
    body.append(arrow(xa, RY-GH/2, xb-xa, GH, col, left=lft))
    cx=(xa+xb)/2
    if name: body.append(txt(cx, RY+GH/2+12, name, 7.6, TXT, weight='bold', style='italic'))
    if dom:
        y = RY+GH/2+ (22 if name else 12)
        body.append(txt(cx, y, dom, 6.2, UNIQ if col==UNIQ else SUB))
out.append(g('genes', *body))
# bracket over the E1-unique genes
ua, ub = km(2820.4), km(2823.5)
br=[line(ua, RY-GH/2-12, ub, RY-GH/2-12, UNIQ, 1.0),
    line(ua, RY-GH/2-12, ua, RY-GH/2-5, UNIQ, 1.0),
    line(ub, RY-GH/2-12, ub, RY-GH/2-5, UNIQ, 1.0),
    txt((ua+ub)/2, RY-GH/2-17, 'E1 strain-unique genes', 8.2, UNIQ, weight='bold')]
out.append(g('unique_bracket', *br))
# legend
LX, LY = 486.0, 56.0
leg=[]
for i,(col,lab) in enumerate([(RM,'Type I R-M system'),(UNIQ,'E1 strain-unique (DUF262/DUF1524)'),
                              (WYL,'WYL phage-defense regulator'),(TA,'Toxin–antitoxin system'),
                              (MOB,'Mobile element / integrase')]):
    y=LY+i*11.5
    leg.append(rect(LX, y-5.6, 8.0, 7.6, fill=col))
    leg.append(txt(LX+11.5, y, lab, 5.9, TXT, anchor='start'))
out.append(g('legend', *leg))
# axis
AXY=132.0
ax=[line(X0-4, AXY, X1+4, AXY, TXT, 0.9)]
for kb,lab in [(2815,'2815 kb'),(2830,'2830 kb')]:
    x=km(kb); ax+=[line(x,AXY,x,AXY+3.4,TXT,0.9), txt(x,AXY+12,lab,6.6,TXT)]
ax.append(txt((X0+X1)/2, AXY+27, 'E1 genomic position (bp)', 8.0, TXT))
out.append(g('axis', *ax))
out.append(tail())
open('svgout/Fig. S5.svg','w',encoding='utf-8').write(''.join(out))
print('S5 ok')
