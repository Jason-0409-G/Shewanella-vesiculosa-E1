#!/usr/bin/env python3
"""Fig. 1b - crop the TEM field shown in the panel from the raw micrograph.

The panel is a 1:1 crop of the preview image, without resampling, rescaling or
level adjustment. The scale bar and labels were added afterwards in the vector
editor.

Input   tem_011_preview.png   (01_raw_data/TEM_microscopy/)
Output  Fig1b_TEM.png
"""
from pathlib import Path
from PIL import Image

SRC = Path("tem_011_preview.png")
OUT = Path("Fig1b_TEM.png")
BOX = (309, 21, 309 + 712, 21 + 1053)   # left, upper, right, lower

im = Image.open(SRC)
im.crop(BOX).save(OUT)
print(f"{SRC} {im.size} -> {OUT} {(BOX[2]-BOX[0], BOX[3]-BOX[1])}")
