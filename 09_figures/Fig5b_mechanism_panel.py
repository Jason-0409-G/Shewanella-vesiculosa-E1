#!/usr/bin/env python3
"""Fig. 5b - proposed host-microbe pathway of terminal cellobiose processing and its carbon fates.

Draws the complete panel: Antarctic krill with its gut and hepatopancreas routes,
host-derived cellobiose, the E1 beta-glucoside-processing island (outer-membrane
porin, periplasmic GH3_e108, BglT, cytoplasmic GH1 bglB, MFS transporter,
GH3_e227) and the fates of cellulose-derived carbon (assimilation, growth,
respiration, faecal-pellet export). The area at the top left is reserved for
panel a. Everything is vector except the raster illustrations.

Input (not distributed with this repository)
--------------------------------------------
data/fig5img/   raster illustrations (krill, faecal pellets, food web, plankton,
                krill biomass; PNG with alpha). They are third-party artwork and
                are not included; without them the script cannot run.

Output (current directory)
--------------------------
Fig5b_mechanism_panel.{svg,pdf,png}

Arial only; text stays editable in the SVG. The published panel was assembled
with panel a and finished in a vector editor.
"""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle, RegularPolygon, Ellipse, PathPatch
from matplotlib.path import Path as MPath
from PIL import Image
import numpy as np

ICON = Path("data/fig5img")
NAVY = "#182448"; SLATE = "#606C78"; TEAL_D = "#245454"; BOXF = "#EEF2F8"; BOXE = "#486CA8"
CELL_OUT = "#3F6DB5"; CELL_IN = "#8FB0D8"; PERI = "#E7EEF7"; CYTO = "#F5F6EE"; CYTO_E = "#8FB08A"
MEM_F = "#D8EBD5"; MEM_E = "#4E8E4E"; BGLT_F = "#5E9E5A"; BGLT_E = "#2E6630"
GOLD = "#E2AA3C"; GOLD_E = "#906030"; GLC = "#6AAE5E"; GLC_E = "#3E7A3E"; RED = "#C8285A"
PATH_B = "#8FB0D8"; AMBER = "#D89A20"; GREY = "#C4C4C4"; MFS_T = "#6A7F99"; GUT = "#4A78C0"
mpl.rcParams.update({"font.family": ["Arial"], "font.size": 12, "pdf.fonttype": 42, "svg.fonttype": "none",
                     "text.color": NAVY, "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                     "mathtext.bf": "Arial:bold"})
DPI = 300
FW, H = 6.4, 5.9                    # fates block (inches)
CW = 5.87                           # cell block (inches); cell units: 150 per inch, y downwards
W, HT = 17.35, 11.3                # whole panel B canvas; the lower block keeps y in [0, H]
fig = plt.figure(figsize=(W, HT), dpi=DPI)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, HT); ax.axis("off")
DASH = (0, (3, 2)); FDASH = (0, (2.2, 1.6))

def F(x, y):                        # fates coords (0-100, 7.6-100) -> inches
    return x / 100 * FW, (y - 7.6) / 92.4 * H
def C(x, y):                        # cell coords (screen-like units, y down) -> inches
    return FW + (x - 900) / 150, (882 - y) / 150

# ---------------- generic helpers (inches) ----------------
def icon(name, p, width_in):
    im = np.asarray(Image.open(ICON / name).convert("RGBA"))
    ax.add_artist(AnnotationBbox(OffsetImage(im, zoom=width_in * 72 / im.shape[1]), p, frameon=False, zorder=3))

def arrow(p0, p1, color=TEAL_D, lw=2.4, ls="-", ms=16, z=2):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, color=color, lw=lw,
                                 ls=ls, shrinkA=0, shrinkB=0, zorder=z, joinstyle="miter"))

def poly(pts, color, lw, ls, z=2):
    xs, ys = zip(*pts); ax.plot(xs, ys, color=color, lw=lw, ls=ls, zorder=z, solid_joinstyle="miter")

def route(pts, color, lw=1.8, ls=DASH, ms=14, z=2, head=True):
    """Polyline through pts (inches), with an arrow head on the last segment."""
    if len(pts) > 2: poly(pts[:-1], color, lw, ls, z)
    if head: arrow(pts[-2], pts[-1], color=color, lw=lw, ls=ls, ms=ms, z=z)
    else: poly(pts[-2:], color, lw, ls, z)

def text(p, s, size=12, ha="center", color=NAVY, weight="bold", **kw):
    ax.text(*p, s, fontsize=size, ha=ha, va="center", color=color, fontweight=weight, zorder=6, **kw)

def rrect_path(x0, y0, x1, y1, r):  # cell units, per-corner radii r=(tl,tr,br,bl); y down
    tl, tr, br, bl = r; k = 0.5523
    v = [(x0 + tl, y0), (x1 - tr, y0)]; c = [MPath.MOVETO, MPath.LINETO]
    def corner(cx, cy, p1, p2):
        v.extend([p1, p2, (cx, cy)]); c.extend([MPath.CURVE4] * 3)
    if tr: corner(x1, y0 + tr, (x1 - tr + k * tr, y0), (x1, y0 + tr - k * tr))
    v.append((x1, y1 - br)); c.append(MPath.LINETO)
    if br: corner(x1 - br, y1, (x1, y1 - br + k * br), (x1 - br + k * br, y1))
    v.append((x0 + bl, y1)); c.append(MPath.LINETO)
    if bl: corner(x0, y1 - bl, (x0 + bl - k * bl, y1), (x0, y1 - bl + k * bl))
    v.append((x0, y0 + tl)); c.append(MPath.LINETO)
    if tl: corner(x0 + tl, y0, (x0, y0 + tl - k * tl), (x0 + tl - k * tl, y0))
    v.append(v[0]); c.append(MPath.CLOSEPOLY)
    return MPath([C(*p) for p in v], c)

def cbox(cx, cy, w, h, s, fc, ec, ls="-", lw=2.0, tcolor=NAVY, size=12, rad=10):
    x0, y1 = C(cx - w / 2, cy + h / 2)
    ax.add_patch(FancyBboxPatch((x0, y1), w / 150, h / 150, boxstyle=f"round,pad=0,rounding_size={rad / 150}",
                                fc=fc, ec=ec, lw=lw, ls=ls, zorder=4))
    if s: text(C(cx, cy), s, size=size, color=tcolor)

def hexa(cx, cy, face, edge, r=24):
    ax.add_patch(RegularPolygon(C(cx, cy), 6, radius=r / 150, orientation=0, fc=face, ec=edge, lw=1.6, zorder=5))

def barrel(cx, cy, w, h, fc, ec, alpha=1.0):
    cbox(cx, cy, w, h, "", fc, ec, lw=2.2, rad=24)
    ax.patches[-1].set_alpha(alpha); ax.patches[-1].set_zorder(5)
    for dx in (-w / 5, 0, w / 5):
        ax.plot(*zip(C(cx + dx, cy - h / 2 + 4), C(cx + dx, cy + h / 2 - 4)), color=ec, lw=1.2, zorder=5, alpha=alpha)

# LEFT: fates of cellulose-derived carbon
text(F(46, 96.5), "Fates of cellulose-derived carbon", size=15)

# assimilation box
x0, y0 = F(45, 76); x1, y1 = F(71, 83.6)
ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0.026,rounding_size=0.1",
                            fc=BOXF, ec=BOXE, lw=1.5, zorder=3))
text(F(58, 81.7), "Assimilation by"); text(F(58, 77.9), "krill holobiont")

# left column: egested cellulose -> faecal pellet -> export
icon("icon_05.png", F(12, 70.5), 1.35)
route([F(8, 62.5), F(8, 55)], SLATE)
text(F(10.5, 60.0), "within krill", 10.5, "left"); text(F(10.5, 56.8), "faecal pellet", 10.5, "left")
icon("icon_01.png", F(12, 46.5), 1.35)
route([F(8, 38.5), F(8, 31)], SLATE)
text(F(10.5, 36.0), "fast-sinking", 10.5, "left"); text(F(10.5, 32.8), "faecal pellet", 10.5, "left")
bx, by = F(2, 12); bx1, by1 = F(22, 29)
ax.add_patch(FancyBboxPatch((bx, by), bx1 - bx, by1 - by, boxstyle="round,pad=0,rounding_size=0.08", fc="#8FB6CF", ec="none", zorder=2))
ax.add_patch(Rectangle((bx, by), bx1 - bx, F(0, 17)[1] - by, fc="#6B5A47", ec="none", zorder=2))
ax.add_patch(Rectangle((bx, by), bx1 - bx, F(0, 13.3)[1] - by, fc="#2E4466", ec="none", zorder=2))
for yy, s in ((25.6, "deep-ocean"), (22.9, "export / seafloor"), (20.2, "sedimentation")):
    text(F(12, yy), s, 8.5, color="white")

# growth: -> krill biomass -> food web
arrow(F(53, 75.3), F(53, 64))
text(F(51.5, 70), "growth", 11, "right", color=TEAL_D)
icon("icon_04.png", F(52, 57), 1.9)
text(F(51, 48.3), "krill biomass", 11)
arrow(F(52, 45.3), F(50, 37.5))
icon("icon_02.png", F(49, 31.5), 1.4)
icon("icon_03.png", F(49, 21), 1.75)
text(F(49, 12), "food web", 12.5)

# respiration: CO2 -> DIC -> phytoplankton -> re-fixation -> krill
arrow(F(64, 75.3), F(77.5, 66.5))
text(F(73.5, 73.8), "respiration", 11, "left", color=TEAL_D)
text(F(84, 62.5), r"CO$_\mathbf{2}$", 13)
route([F(84, 59), F(84, 54)], SLATE)
text(F(84, 50.8), "DIC pool", 12.5)
route([F(84, 47.5), F(84, 42)], SLATE)
icon("icon_06.png", F(86, 28.5), 1.15)
route([F(77, 26), F(71, 26), F(71, 57), F(67.8, 57)], SLATE)
text(F(69.5, 40.5), "photosynthetic", 9.5, "right", color=SLATE, weight="normal")
text(F(69.5, 37.3), "re-fixation", 9.5, "right", color=SLATE, weight="normal")

# RIGHT: E1 beta-glucoside-processing island
text(C(1310, 30), "E1 β-glucoside-processing island", size=15)
# envelope: outer membrane (thick) / white gap / inner line; periplasm fill
ax.add_patch(PathPatch(rrect_path(940, 108, 1690, 850, (110,) * 4), fc="white", ec=CELL_OUT, lw=6, zorder=1))
ax.add_patch(PathPatch(rrect_path(954, 122, 1676, 836, (96,) * 4), fc=PERI, ec=CELL_IN, lw=2.2, zorder=1))
# cytoplasmic compartment + inner membrane
ax.add_patch(PathPatch(rrect_path(962, 130, 1252, 828, (88, 0, 140, 88)), fc=CYTO, ec=CYTO_E, lw=1.6, zorder=1))
mx0, my0 = C(1250, 828); mx1, my1 = C(1312, 130)
ax.add_patch(Rectangle((mx0, my0), mx1 - mx0, my1 - my0, fc=MEM_F, ec=MEM_E, lw=2.2, zorder=2))
for yy in np.arange(142, 822, 13):
    if 218 < yy < 290 or 618 < yy < 692: continue
    for xx, tail in ((1265, 1277), (1297, 1285)):
        ax.plot(*zip(C(xx, yy), C(tail, yy)), color=MEM_E, lw=0.5, zorder=2)
        ax.add_patch(plt.Circle(C(xx, yy), 4.2 / 150, fc=MEM_F, ec=MEM_E, lw=0.6, zorder=3))

# transporters
barrel(1280, 253, 110, 60, BGLT_F, BGLT_E)
text(C(1322, 205), "BglT", 14, "left")
barrel(1275, 655, 106, 62, "#D6D6D6", "#AFAFAF", alpha=0.85)
text(C(1360, 705), "MFS", 14, "center", color=MFS_T)
# porin in the outer membrane + host cellobiose entering from the right
px, py = C(1664, 560)
ax.add_patch(FancyBboxPatch((px, py), 44 / 150, 140 / 150, boxstyle=f"round,pad=0,rounding_size={20 / 150}",
                            fc="#5B8CC8", ec="#2F5E9E", lw=2.2, zorder=5))
ax.add_patch(Ellipse(C(1686, 440), 28 / 150, 13 / 150, fc="#A8C4E6", ec="#2F5E9E", lw=1.4, zorder=6))

# cellobiose (periplasm and cytoplasm), glucose
for cx in (1137, 1463):
    hexa(cx, 228, GOLD, GOLD_E); hexa(cx, 276, GOLD, GOLD_E)
hexa(1378, 655, GLC, GLC_E); hexa(1137, 662, GLC, GLC_E)

# enzymes
cbox(1137, 467, 150, 54, "GH1 bglB", "#EEF5EC", "#5E9E5A", lw=2.2)
cbox(1483, 453, 150, 56, "GH3_e108", "#F2F2F2", GREY, ls=(0, (3, 2)), lw=1.8)
cbox(1497, 652, 162, 58, "GH3_e227", "#FFF3C4", "#E0B040", lw=2.2)
cbox(1000, 385, 62, 42, r"CO$_\mathbf{2}$", "#EEF7EC", "#6AAE5E", lw=1.6, size=12, rad=8)
text(C(1026, 548), "central carbon\nmetabolism", 11, rotation=-90, linespacing=1.05)

# pathways
route([C(1612, 500), C(1612, 253), C(1492, 253)], RED, lw=1.8, ls=DASH, head=False)       # porin -> periplasmic cellobiose
route([C(1437, 253), C(1337, 253)], RED, lw=2.0, ls=DASH)                              # -> BglT
route([C(1664, 500), C(1612, 500), C(1612, 655), C(1580, 655)], AMBER, lw=1.8, ls=FDASH, head=False)
route([C(1560, 453), C(1612, 453)], GREY, lw=1.6, ls=FDASH, head=False)
route([C(1408, 453), C(1378, 453), C(1378, 630)], GREY, lw=1.6, ls=FDASH, head=False)
route([C(1353, 655), C(1330, 655)], PATH_B, lw=1.8)                                    # glucose -> MFS
route([C(1221, 662), C(1164, 662)], PATH_B, lw=1.8)                                    # MFS -> cytoplasm
route([C(1224, 253), C(1164, 253)], PATH_B, lw=1.8)                                    # BglT -> cytoplasm
route([C(1137, 304), C(1137, 438)], PATH_B, lw=1.8)
route([C(1137, 496), C(1137, 636)], PATH_B, lw=1.8)
route([C(1112, 662), C(1000, 662), C(1000, 408)], PATH_B, lw=1.8)                      # -> central C metabolism -> CO2

# cross-links
GUT_END = C(1652, 118)             # the gut route arrives at the top-right corner of the island
text(F(10.5, 87.6), "unhydrolysed", 10, "left", color=SLATE, weight="normal")
text(F(10.5, 85.2), "cellulose", 10, "left", color=SLATE, weight="normal")

# TOP: Antarctic krill and the two digestive routes
text((10.75, 10.95), "B", 22, "left", color="#231F20")
text((14.0, 10.62), "Antarctic krill", 20)
text((14.0, 10.08), "(Euphausia superba)", 13, style="italic")
KX0, KW, KTOP = 10.75, 6.55, 9.75
kim = np.asarray(Image.open(ICON / "krill_big_0.png").convert("RGBA")); KH = KW * kim.shape[0] / kim.shape[1]
ax.imshow(kim, extent=(KX0, KX0 + KW, KTOP - KH, KTOP), zorder=3, interpolation="lanczos")
ax.set_xlim(0, W); ax.set_ylim(0, HT)
gut0 = (KX0 + 0.395 * KW, KTOP - 0.15 * KH)            # mid-gut
hep0 = (KX0 + 0.615 * KW, KTOP - 0.37 * KH)            # hepatopancreas
GD = (0, (4, 2.5))
route([gut0, GUT_END], GUT, lw=2.0, ls=GD, ms=18, z=4)
text(((gut0[0] + GUT_END[0]) / 2 + 0.45, (gut0[1] + GUT_END[1]) / 2 + 0.2), "Gut", 14, "left")
route([hep0, (hep0[0] - 0.08, 6.08)], GUT, lw=2.0, ls=GD, ms=18, z=4)
text((hep0[0] + 0.25, 7.45), "Hepatopancreas", 14, "left")
# unhydrolysed cellulose branches off the gut route and is egested without assimilation
yb = F(0, 90.5)[1]; tb = (gut0[1] - yb) / (gut0[1] - GUT_END[1]); xb = gut0[0] + tb * (GUT_END[0] - gut0[0])
route([(xb, yb), F(8, 90.5), F(8, 80.5)], SLATE)

# RIGHT: host-derived cellobiose
HX0, CWD, CG = 12.35, 0.66, 0.018
hx = [HX0 + i * (CWD + CG) for i in range(6)]
text(((hx[0] + hx[-1] + CWD) / 2, 5.68), "Host-derived cellobiose", 15)
EPI_F, EPI_E, NUC = "#EAF1F9", "#A8C4E6", "#7FA3D0"
def epi_row(y0, villi_top):
    for x in hx:
        ax.add_patch(FancyBboxPatch((x, y0), CWD, 0.6, boxstyle="round,pad=0,rounding_size=0.09",
                                    fc=EPI_F, ec=EPI_E, lw=1.6, zorder=3))
        ax.add_patch(Ellipse((x + CWD / 2, y0 + 0.3), 0.13, 0.09, fc=NUC, ec="none", zorder=4))
        yv = y0 + 0.6 if villi_top else y0
        for j in range(4):
            cx = x + 0.17 + j * 0.105
            th = np.linspace(0, np.pi, 30); sgn = 1 if villi_top else -1
            ax.plot(cx + 0.05 * np.cos(th), yv + sgn * 0.07 * np.sin(th), color=EPI_E, lw=1.3, zorder=3)
epi_row(4.62, villi_top=False)
epi_row(0.28, villi_top=True)
# cellulose chain (host lumen)
CX = hx[-1] + CWD / 2; CEL_F, CEL_E = "#62AFA8", "#1E5E58"; RH = 0.16
ys = [4.08 - i * 0.39 for i in range(5)]
for y in ys:
    ax.add_patch(RegularPolygon((CX, y), 6, radius=RH, orientation=0, fc=CEL_F, ec=CEL_E, lw=2.0, zorder=5))
for y0_, y1_ in zip(ys[:-1], ys[1:]):
    ax.plot([CX, CX], [y0_ - RH, y1_ + RH], color=CEL_E, lw=3, zorder=4)
for k in range(3):
    ax.add_patch(plt.Circle((CX, ys[-1] - RH - 0.11 - k * 0.12), 0.035, color=CEL_E, zorder=5))
text((CX + 0.42, 3.3), "cellulose", 15, rotation=-90)
# cellobiose
BX, BY = 12.66, 2.69
for dy in (0.2, -0.2):
    ax.add_patch(RegularPolygon((BX, BY + dy), 6, radius=0.2, orientation=0, fc="#E6AC3C", ec="#9A6424", lw=2.0, zorder=5))
ax.plot([BX, BX], [BY - 0.02, BY + 0.02], color="#9A6424", lw=4, zorder=4)
text((BX + 0.4, BY), "cellobiose", 15, rotation=-90)
# host enzymes: cellulose -> cellobiose
arrow((CX - 0.3, BY), (BX + 0.62, BY), color="#285858", lw=2.2, ms=18)
EX = (CX - 0.3 + BX + 0.62) / 2
text((EX, 3.36), "host", 13)
text((EX, 2.97), "endo-β-1,4-glucanase", 13, style="italic")
text((EX, 2.42), "exo-β-1,4-glucanase", 13, style="italic")
# cellobiose -> outer-membrane porin of E1
arrow((BX - 0.21, BY), C(1712, 478), color="#285858", lw=2.4, ms=16)

# glucose released by E1 hydrolysis -> holobiont assimilation
gy = F(0, 79.8)[1]
route([(C(938, 0)[0], gy), (F(71.8, 0)[0], gy)], SLATE)
text(((C(938, 0)[0] + F(71.8, 0)[0]) / 2, gy + 0.16), "glucose", 10, color=SLATE, weight="normal")

for ext in ("svg", "pdf", "png"):
    fig.savefig(f"Fig5b_mechanism_panel.{ext}", dpi=DPI)
print("written: Fig5b_mechanism_panel.{svg,pdf,png}")
