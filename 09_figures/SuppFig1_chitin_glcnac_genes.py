"""Fig. S1 - chitin and GlcNAc utilisation genes of S. vesiculosa E1.

Hand-laid arrow map in editable SVG. Left: gbpA and chiA (E1 chromosome,
0.1-1.0 Mb); right: the contiguous nagA-glmS-glucosamine kinase-chb cluster at
3.888-3.895 Mb. Gene coordinates are the E1 Prokka annotation (L5.gff).
Run from 09_figures/; writes svgout/Fig. S1.svg.
"""
import os, sys; sys.path.insert(0, '.')
os.makedirs('svgout', exist_ok=True)
from _svglib import *

H=150.0
BLUE='#1F4E79'; TXT='#333333'; GREY='#8A8A8A'; RED='#C62828'; PINK='#FCF3F3'; RAIL='#E0E0E0'
RY=56.0; GH=20.0
NY, SY, CY = RY+GH/2+11.5, RY+GH/2+21.0, RY+GH/2+30.0   # text rows: name / annotation / coordinate
AXY=122.0
out=[head(H,'Fig. S1 Chitin and GlcNAc utilisation genes')]

# left half
LX0, LX1 = 20.0, 222.0
def lmap(mb): return LX0 + (mb-0.10)/(1.00-0.10)*(LX1-LX0)
left=[line(LX0, RY, LX1, RY, RAIL, 3.4, cap='round')]
for name, sub, mb, wid in [('gbpA','AA10 LPMO + CBM73',0.160,44.0),
                           ('chiA','GH18 + CBM5 (endo-chitinase)',0.889,56.0)]:
    cx=lmap(mb); left.append(arrow(cx-wid/2, RY-GH/2, wid, GH, BLUE))
    left.append(txt(cx, NY, name, 8, TXT, weight='bold', style='italic'))
    left.append(txt(cx, SY, sub, 5.2, '#555', style='italic'))
    left.append(txt(cx, CY, f'{mb:.3f} Mb', 5.2, GREY))
left.append(line(LX0, AXY, LX1, AXY, TXT, 0.9))
for mb,lab in [(0.2,'0.2 Mb'),(0.5,'0.5 Mb'),(0.8,'0.8 Mb')]:
    x=lmap(mb); left+=[line(x,AXY,x,AXY+3.2,TXT,0.9), txt(x,AXY+11,lab,6.2,TXT)]
out.append(g('left_region', *left))

# axis break
BX=248.0
brk=[line(BX+dx-3.4, RY+8, BX+dx+3.4, RY-8, GREY, 1.3) for dx in (-6.5,-2,2,6.5)]
brk.append(txt(BX, AXY+11, '~3 Mb gap', 6.2, GREY, style='italic'))
out.append(g('gap_break', *brk))

# right half
BL, BR = 272.0, 583.0
out_box=rect(BL, 10.0, BR-BL, 130.0, fill=PINK, stroke=RED, sw=1.0, rx=7, dash='4 3')
RX0, RX1 = BL+18, BR-14
K0, K1 = 3888.5, 3895.5
def kmap(kb): return RX0 + (kb-K0)/(K1-K0)*(RX1-RX0)
# (name, annotation, coordinate, start kb, end kb, name font size)
genes=[('nagA','GlcNAc-6-P deacetylase','3,888,921–3,890,072',3888.921,3890.072, 7.4),
       ('glmS','GlcN-6-P synthase','3,890,109–3,891,110',3890.109,3891.110, 7.4),
       ('glucosamine kinase','EC 2.7.1.8 (ROK/BadF family)','3,891,193–3,892,113',3891.193,3892.113, 6.0),
       ('chb','GH20 exo-hexosaminidase','3,892,371–3,895,061',3892.371,3895.061, 7.4)]
right=[out_box,
       txt((BL+BR)/2, 26.0, 'contiguous co-directional cluster (37–258 bp gaps)', 7.2, RED, weight='bold'),
       line(RX0, RY, RX1, RY, RAIL, 3.4, cap='round')]
for name, sub, coord, a, b, fs in genes:
    xa, xb = kmap(a), kmap(b); cx=(xa+xb)/2
    right.append(arrow(xa, RY-GH/2, xb-xa, GH, BLUE))
    right.append(txt(cx, NY, name, fs, TXT, weight='bold', style='italic'))
    right.append(txt(cx, SY, sub, 4.4, '#555', style='italic'))
    right.append(txt(cx, CY, coord, 4.4, GREY))
right.append(line(RX0, AXY, RX1, AXY, TXT, 0.9))
for kb,lab in [(3890,'3890 kb'),(3895,'3895 kb')]:
    x=kmap(kb); right+=[line(x,AXY,x,AXY+3.2,TXT,0.9), txt(x,AXY+11,lab,6.2,GREY)]
out.append(g('cluster_region', *right))
out.append(tail())
open('svgout/Fig. S1.svg','w',encoding='utf-8').write(''.join(out))
print('S1 ok')
