"""Fig. S2 - CAZyme repertoire of the eight Shewanella genomes (Sankey).

Per-strain counts of the six CAZy classes (GH, GT, PL, CE, AA, CBM) from the
dbCAN3 >=2-of-3-tool consensus; the numbers live in SuppFig2_cazyme_sankey_data.py.
E1 ribbons are highlighted and the E1 GH / GT / CBM counts are labelled.
Run from 09_figures/ after `mkdir -p svgout`; writes svgout/Fig. S2.svg.
"""
import sys; sys.path.insert(0,'.')
from _svglib import *
from SuppFig2_cazyme_sankey_data import STRAINS, CLASSES, FLOW

H=430.0
SCOL={'S. vesiculosa E1':'#9F392C','S. vesiculosa M7':'#C4806D','S. frigidimarina':'#6A9FBE',
      'S. sp002836315':'#79AE89','S. sp014164505':'#B0A3C8','S. livingstonensis':'#C4A96B',
      'S. psychromarinicola':'#83B5AE','S. polaris':'#A8A8A8'}
CCOL={'GH':'#48638D','GT':'#699FC6','PL':'#79AE89','CE':'#D9903C','AA':'#9C7BAA','CBM':'#7A7A7A'}
CFULL={'GH':'GH Glycoside Hydrolases','GT':'GT Glycosyl Transferases','PL':'PL Polysaccharide Lyases',
       'CE':'CE Carbohydrate Esterases','AA':'AA Auxiliary Activities','CBM':'CBM Carbohydrate-Binding Modules'}
TXT='#222222'

T, B = 18.0, 412.0
LX, RX = 176.0, 404.0        # left / right x of the node bars
NW = 11.0                    # node bar width
GAPL, GAPR = 6.0, 7.0
TOT=sum(v for _,v in STRAINS)
unit=min((B-T-GAPL*(len(STRAINS)-1))/TOT, (B-T-GAPR*(len(CLASSES)-1))/TOT)

def stack(items, gap):
    tot=sum(v for _,v in items)*unit + gap*(len(items)-1)
    y=T+((B-T)-tot)/2; pos={}
    for k,v in items:
        pos[k]=(y, v*unit); y+=v*unit+gap
    return pos
LP=stack(STRAINS, GAPL); RP=stack(CLASSES, GAPR)

out=[head(H,'Fig. S2 CAZyme repertoire of the eight Shewanella genomes')]

# ribbons (drawn first so the nodes sit on top)
loff={k:LP[k][0] for k in LP}; roff={k:RP[k][0] for k in RP}
ribbons=[]
for sname,_ in STRAINS:
    for cname,_ in CLASSES:
        v=FLOW[f'{sname}|{cname}']
        if not v: continue
        h=v*unit
        y0, y1 = loff[sname], roff[cname]
        x0, x1 = LX+NW, RX
        cx=(x0+x1)/2
        d=(f'M {x0:.2f} {y0:.2f} C {cx:.2f} {y0:.2f} {cx:.2f} {y1:.2f} {x1:.2f} {y1:.2f} '
           f'L {x1:.2f} {y1+h:.2f} C {cx:.2f} {y1+h:.2f} {cx:.2f} {y0+h:.2f} {x0:.2f} {y0+h:.2f} Z')
        op = 0.72 if sname=='S. vesiculosa E1' else 0.34
        ribbons.append(path(d, fill=SCOL[sname], op=op))
        loff[sname]+=h; roff[cname]+=h
out.append(g('ribbons', *ribbons))

# strain nodes
ln=[]
for sname, n in STRAINS:
    y,h=LP[sname]
    ln.append(rect(LX, y, NW, h, fill=SCOL[sname]))
    cy=y+h/2
    ln.append(txt(LX-8, cy-1.0, sname, 8.2, TXT, anchor='end', style='italic',
                  weight='bold' if sname=='S. vesiculosa E1' else 'normal'))
    ln.append(txt(LX-8, cy+9.0, f'(n={n})', 7.2, '#555', anchor='end'))
out.append(g('strain_nodes', *ln))

# class nodes
rn=[]
for cname, n in CLASSES:
    y,h=RP[cname]
    rn.append(rect(RX, y, NW, h, fill=CCOL[cname]))
    cy=y+h/2
    if h < 26:      # short node: put name and count on one line so they do not collide
        rn.append(txt(RX+NW+8, cy+2.6, f'{CFULL[cname]}  (n={n})', 7.6, TXT, anchor='start', style='italic'))
    else:
        rn.append(txt(RX+NW+8, cy-1.0, CFULL[cname], 8.2, TXT, anchor='start', style='italic'))
        rn.append(txt(RX+NW+8, cy+9.0, f'(n={n})', 7.2, '#555', anchor='start'))
out.append(g('class_nodes', *rn))

# count labels on the E1 GH / GT / CBM ribbons
lo={k:LP[k][0] for k in LP}; ro={k:RP[k][0] for k in RP}
labels=[]
for sname,_ in STRAINS:
    for cname,_ in CLASSES:
        v=FLOW[f'{sname}|{cname}']
        if not v: continue
        h=v*unit
        if sname=='S. vesiculosa E1' and cname in ('GH','GT','CBM'):
            tt=0.16                     # near the left end, before the ribbons cross
            y0, y1 = lo[sname]+h/2, ro[cname]+h/2
            mx = (LX+NW) + (RX-(LX+NW))*tt
            my = (1-tt)**3*y0 + 3*(1-tt)**2*tt*y0 + 3*(1-tt)*tt**2*y1 + tt**3*y1
            labels.append(txt(mx, my+3.0, str(v), 8.8, '#FFFFFF', weight='bold',
                              stroke=SCOL[sname], sw=2.0))
        lo[sname]+=h; ro[cname]+=h
out.append(g('flow_labels', *labels))
out.append(tail())
open('svgout/Fig. S2.svg','w',encoding='utf-8').write(''.join(out))
print('S2 ok')
