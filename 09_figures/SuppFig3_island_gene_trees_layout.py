"""Fig. S3 - protein ML trees for the four island genes (submitted layout).

Topology, branch lengths and SH-aLRT/UFBoot support values come from the
IQ-TREE .treefile outputs; the newick strings are inlined so the figure
reproduces without those files, and leaf order is fixed to match the
submitted figure. Nodes are labelled when SH-aLRT >= 85.
SuppFig3_island_gene_trees.py is the data-driven counterpart: it reads the
.treefile outputs directly (matplotlib layout, labels at either value >= 80)
and does not reproduce the submitted layout exactly.
Run from 09_figures/ after `mkdir -p svgout`; writes svgout/Fig. S3.svg.
"""
import sys, io, re; sys.path.insert(0,'.')
from _svglib import *
from Bio import Phylo

TREES = {
 'A': ('GH1 (bglB, OG0001451)', 0.1,
  "(L5_PB002_01809:0.0065339073,M7_KDH10_01698:0.0047219197,(((frig_FRIG_01383:0.0258917901,"
  "sp02_SP02_02811:0.0120224756)96/92:0.0220611077,pola_POLA_02651:0.0519148474)0.9/27:0.0071804362,"
  "(livi_LIVI_03054:0.0410873628,psyc_PSYC_04062:0.0550455783)90.1/52:0.0152981797)97.8/95:0.0218676131);"),
 'B': ('GH3_e108 (OG0002997)', 0.1,
  "(L5_PB002_01800:0.0000020975,M7_KDH10_01707:0.0093787896,(((frig_FRIG_01374:0.0056373472,"
  "sp14_SP14_04366:0.0256914840)12.8/47:0.0033242975,sp02_SP02_02820:0.0126264458)100/100:0.1018810724,"
  "(livi_LIVI_03063:0.0360262151,pola_POLA_02660:0.0367043055)99.6/100:0.0458512543)100/100:0.0872068223);"),
 'C': ('GH3_e227 (OG0001452)', 0.05,
  "((S_vesiculosa_E1:0.0060044416,S_vesiculosa_M7:0.0070303400)100/100:0.1055596697,"
  "((((S_frigidimarina:0.0375963427,S_sp002836315:0.0200219122)91.4/90:0.0147721757,"
  "S_sp014164505:0.0219038211)100/100:0.1149432933,S_psychromarinicola:0.0723227712)94.2/89:0.0263443574,"
  "S_polaris:0.0550992552)79.5/76:0.0189034949,S_livingstonensis:0.0641803483);"),
 'D': ('MFS transporter (OG0000506)', 0.05,
  "(PB002_01810:0.0000010000,KDH10_01697:0.0000010000,((((FRIG_01384:0.0230784360,"
  "(SP02_02810:0.0070422514,SP14_03070:0.0075405957)94.4/97:0.0167392640)99.2/92:0.0416665950,"
  "PSYC_04061:0.0268884349)91.8/60:0.0209809706,LIVI_03053:0.0374591989)46/41:0.0096391113,"
  "POLA_02649:0.0504456851)100/100:0.0693355930);"),
}
def disp(name):
    n=name.upper()
    if 'L5' in n or 'PB002' in n or 'VESICULOSA_E1' in n: return 'S. vesiculosa E1'
    if 'M7' in n or 'KDH10' in n or 'VESICULOSA_M7' in n: return 'S. vesiculosa M7'
    if 'FRIG' in n: return 'S. frigidimarina'
    if 'SP02' in n or 'SP002836315' in n: return 'S. sp002836315'
    if 'SP14' in n or 'SP014164505' in n: return 'S. sp014164505'
    if 'LIVI' in n: return 'S. livingstonensis'
    if 'PSYC' in n: return 'S. psychromarinicola'
    if 'POLA' in n: return 'S. polaris'
    return name

E1RED='#A03A2C'; BR='#1A1A1A'; SUP='#7A7A7A'; TXT='#111111'
SUP_MIN=85.0           # SH-aLRT threshold, matching the nodes labelled in the published figure

# leaf order of the published figure (top to bottom); fixed here so the layout matches
LEAF_ORDER={
 'A':['S. livingstonensis','S. psychromarinicola','S. polaris','S. sp002836315','S. frigidimarina',
      'S. vesiculosa M7','S. vesiculosa E1'],
 'B':['S. livingstonensis','S. polaris','S. sp002836315','S. sp014164505','S. frigidimarina',
      'S. vesiculosa M7','S. vesiculosa E1'],
 'C':['S. livingstonensis','S. polaris','S. psychromarinicola','S. sp014164505','S. sp002836315',
      'S. frigidimarina','S. vesiculosa M7','S. vesiculosa E1'],
 'D':['S. livingstonensis','S. psychromarinicola','S. sp014164505','S. sp002836315','S. frigidimarina',
      'S. polaris','S. vesiculosa M7','S. vesiculosa E1'],
}
def build(nwk, key):
    t=Phylo.read(io.StringIO(nwk),'newick')
    og=next(c for c in t.get_terminals() if disp(c.name)=='S. livingstonensis')
    t.root_with_outgroup(og)
    rank={n:i for i,n in enumerate(LEAF_ORDER[key])}
    def keyf(c):
        return min(rank[disp(l.name)] for l in c.get_terminals())
    def sort(c):
        if c.is_terminal(): return
        c.clades.sort(key=keyf)
        for k in c.clades: sort(k)
    sort(t.root)
    return t

def draw_panel(key, ox, oy, pw, ph):
    title, scale, nwk = TREES[key]
    t=build(nwk, key)
    leaves=t.get_terminals()
    depths=t.depths(unit_branch_lengths=False)
    maxd=max(depths[l] for l in leaves) or 1.0
    LAB_W=104.0
    x0=ox+13.0; x1=ox+pw-LAB_W
    ytop=oy+30.0; ybot=oy+ph-30.0
    dy=(ybot-ytop)/max(len(leaves)-1,1)
    ypos={}
    for i,l in enumerate(leaves): ypos[l]=ytop+i*dy
    def yof(c):
        if c.is_terminal(): return ypos[c]
        return sum(yof(k) for k in c.clades)/len(c.clades)
    def xof(c): return x0 + depths[c]/maxd*(x1-x0)
    segs=[]; labs=[]; sups=[]
    def walk(c):
        cy=yof(c); cx=xof(c)
        if c.is_terminal():
            nm=disp(c.name); red = nm=='S. vesiculosa E1'
            labs.append(txt(cx+4.0, cy+2.6, nm, 7.6, E1RED if red else TXT,
                            anchor='start', style='italic', weight='bold' if red else 'normal'))
            return
        ys=[yof(k) for k in c.clades]
        segs.append(line(cx, min(ys), cx, max(ys), BR, 0.9))      # vertical connector
        for k in c.clades:
            segs.append(line(cx, yof(k), xof(k), yof(k), BR, 0.9))  # horizontal branch
            walk(k)
        cf=str(c.confidence) if c.confidence is not None else (c.name or '')
        m=re.match(r'^([\d.]+)/([\d.]+)$', cf or '')
        if m and float(m.group(1))>=SUP_MIN:
            sups.append(txt(cx-2.0, cy+8.4, cf, 5.8, SUP, anchor='end'))
    root=t.root
    segs.append(line(x0-6.0, yof(root), x0, yof(root), BR, 0.9))
    walk(root)
    sx0=ox+13.0; sw=scale/maxd*(x1-x0); sy=oy+ph-12.0
    sc=[line(sx0, sy, sx0+sw, sy, BR, 1.2), txt(sx0+sw/2, sy+10.0, f'{scale}', 7.0, TXT)]
    head_=[txt(ox, oy+11.0, key, 11.0, TXT, anchor='start', weight='bold'),
           txt(ox+13.0, oy+24.0, title, 8.6, TXT, anchor='start', weight='bold')]
    return g(f'panel_{key}', *head_, *segs, *sups, *labs, *sc)

W_, H_ = 595.276, 580.0
PW, PH = 292.0, 278.0
out=[head(H_,'Fig. S3 ML protein phylogenies of the four island genes')]
for key,(ox,oy) in zip('ABCD', [(8,6),(300,6),(8,292),(300,292)]):
    out.append(draw_panel(key, ox, oy, PW, PH))
out.append(tail())
open('svgout/Fig. S3.svg','w',encoding='utf-8').write(''.join(out))
print('S3 ok')
