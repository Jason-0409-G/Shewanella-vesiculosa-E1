# -*- coding: utf-8 -*-
"""Minimal helpers for editable SVG output: text is always emitted as <text>,
shapes are grouped."""
W = 595.276          # 210 mm, A4 width

def head(h, title):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{W}pt" height="{h}pt" viewBox="0 0 {W} {h}" version="1.1">\n'
            f'<title>{title}</title>\n'
            f'<style>text{{font-family:Arial,Helvetica,sans-serif;white-space:pre}}</style>\n')

def tail(): return '</svg>\n'

def g(id_, *body, transform=None):
    t=f' transform="{transform}"' if transform else ''
    return f'<g id="{id_}"{t}>\n' + ''.join(body) + '</g>\n'

def txt(x, y, s, size=7, fill='#333', anchor='middle', weight='normal', style='normal', ls=None,
        stroke=None, sw=0.0):
    e=s.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')
    a=f' letter-spacing="{ls}"' if ls else ''
    k=(f' stroke="{stroke}" stroke-width="{sw}" paint-order="stroke fill" '
       f'stroke-linejoin="round"') if stroke else ''
    return (f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" '
            f'font-size="{size}" fill="{fill}" '
            f'text-anchor="{anchor}" font-weight="{weight}" font-style="{style}"{a}{k}>{e}</text>\n')

def line(x1,y1,x2,y2,stroke='#333',w=0.8,dash=None,cap='butt'):
    d=f' stroke-dasharray="{dash}"' if dash else ''
    return (f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{stroke}" stroke-width="{w}" stroke-linecap="{cap}"{d}/>\n')

def rect(x,y,w,h,fill='#fff',stroke='none',sw=0.8,rx=0,dash=None):
    d=f' stroke-dasharray="{dash}"' if dash else ''
    return (f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>\n')

def path(d, fill='none', stroke='none', sw=0.8, op=1.0):
    return f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" fill-opacity="{op}"/>\n'

def arrow(x, y, w, h, fill, left=True, tip=None):
    """Gene arrow: x, y is the top-left corner; left=True points the arrow leftwards."""
    tip = tip if tip is not None else min(h*0.62, w*0.45)
    if left:
        d=(f'M {x+tip:.2f} {y:.2f} L {x+w:.2f} {y:.2f} L {x+w:.2f} {y+h:.2f} '
           f'L {x+tip:.2f} {y+h:.2f} L {x:.2f} {y+h/2:.2f} Z')
    else:
        d=(f'M {x:.2f} {y:.2f} L {x+w-tip:.2f} {y:.2f} L {x+w:.2f} {y+h/2:.2f} '
           f'L {x+w-tip:.2f} {y+h:.2f} L {x:.2f} {y+h:.2f} Z')
    return path(d, fill=fill)
