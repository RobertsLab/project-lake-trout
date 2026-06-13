#!/usr/bin/env python3
"""Build the lake trout PAV + methylation ecotype deck."""
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION, XL_LABEL_POSITION

# ---------- palette (Ocean Gradient: deep water -> shallow) ----------
DEEP    = RGBColor(0x06, 0x2A, 0x42)   # near-black deep water (title bg)
MID     = RGBColor(0x06, 0x5A, 0x82)   # deep blue
TEAL    = RGBColor(0x1C, 0x72, 0x93)   # teal
SEAFOAM = RGBColor(0x2A, 0x9D, 0x8F)   # green accent (methylation)
CORAL   = RGBColor(0xE9, 0x96, 0x4E)   # warm accent (PAV / contrast)
ICE     = RGBColor(0xCA, 0xDC, 0xFC)   # ice blue light
WHITE   = RGBColor(0xFF, 0xFF, 0xFF)
INK     = RGBColor(0x1A, 0x2B, 0x35)   # body text on light
MUTE    = RGBColor(0x5B, 0x70, 0x7B)   # muted caption
LIGHTBG = RGBColor(0xF4, 0xF7, 0xF9)   # light content bg

HEAD_FONT = "Georgia"
BODY_FONT = "Calibri"

EMU_W, EMU_H = Inches(13.333), Inches(7.5)

prs = Presentation()
prs.slide_width = EMU_W
prs.slide_height = EMU_H
BLANK = prs.slide_layouts[6]


def slide(bg=WHITE):
    s = prs.slides.add_slide(BLANK)
    fill = s.background.fill
    fill.solid()
    fill.fore_color.rgb = bg
    return s


def rect(s, x, y, w, h, color, line=None, shape=MSO_SHAPE.RECTANGLE, shadow=False):
    sp = s.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    sp.fill.solid()
    sp.fill.fore_color.rgb = color
    if line is None:
        sp.line.fill.background()
    else:
        sp.line.color.rgb = line
        sp.line.width = Pt(1)
    sp.shadow.inherit = False
    if shadow:
        el = sp._element.spPr
        # lightweight soft shadow via XML
        from pptx.oxml.ns import qn
        efflst = el.makeelement(qn('a:effectLst'), {})
        outer = el.makeelement(qn('a:outerShdw'),
                               {'blurRad': '60000', 'dist': '25000', 'dir': '5400000', 'rotWithShape': '0'})
        clr = el.makeelement(qn('a:srgbClr'), {'val': '0A2230'})
        alpha = el.makeelement(qn('a:alpha'), {'val': '28000'})
        clr.append(alpha)
        outer.append(clr)
        efflst.append(outer)
        el.append(efflst)
    return sp


def txt(s, x, y, w, h, runs, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP,
        wrap=True, space_after=None):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = wrap
    tf.vertical_anchor = anchor
    tf.margin_left = 0
    tf.margin_right = 0
    tf.margin_top = 0
    tf.margin_bottom = 0
    if isinstance(runs, str):
        runs = [[(runs, {})]]
    for i, para in enumerate(runs):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        if space_after is not None:
            p.space_after = Pt(space_after)
        if isinstance(para, tuple):
            para = [para]
        for text, opt in para:
            r = p.add_run()
            r.text = text
            f = r.font
            f.name = opt.get("font", BODY_FONT)
            f.size = Pt(opt.get("size", 16))
            f.bold = opt.get("bold", False)
            f.italic = opt.get("italic", False)
            f.color.rgb = opt.get("color", INK)
            if opt.get("spacing"):
                _set_spacing(r, opt["spacing"])
    return tb


def _set_spacing(run, pts):
    from pptx.oxml.ns import qn
    rPr = run._r.get_or_add_rPr()
    rPr.set('spc', str(int(pts * 100)))


def bullets(s, x, y, w, h, items, size=16, color=INK, gap=10, bullet_color=TEAL):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = 0
    tf.margin_right = 0
    tf.margin_top = 0
    for i, item in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(gap)
        lead = "—  "
        r = p.add_run(); r.text = lead
        r.font.name = BODY_FONT; r.font.size = Pt(size); r.font.bold = True
        r.font.color.rgb = bullet_color
        if isinstance(item, str):
            item = [(item, {})]
        for text, opt in item:
            r = p.add_run(); r.text = text
            r.font.name = BODY_FONT
            r.font.size = Pt(opt.get("size", size))
            r.font.bold = opt.get("bold", False)
            r.font.italic = opt.get("italic", False)
            r.font.color.rgb = opt.get("color", color)
    return tb


def kicker(s, x, y, text, color=CORAL):
    txt(s, x, y, 8, 0.35, [[(text.upper(), {"font": BODY_FONT, "size": 13, "bold": True,
                                            "color": color, "spacing": 2.5})]])


def page_no(s, n):
    txt(s, 12.4, 7.05, 0.7, 0.3, [[(f"{n:02d}", {"size": 11, "color": MUTE})]],
        align=PP_ALIGN.RIGHT)


N = 0
def num():
    global N
    N += 1
    return N

# ============================================================ 1. TITLE
s = slide(DEEP)
# layered depth bands
rect(s, 0, 5.0, 13.333, 2.5, MID)
rect(s, 0, 6.0, 13.333, 1.5, TEAL)
rect(s, 0, 6.8, 13.333, 0.7, SEAFOAM)
kicker(s, 0.9, 1.05, "Salvelinus namaycush  ·  Roberts Lab, UW SAFS", ICE)
txt(s, 0.9, 1.6, 11.5, 2.4, [
    [("Two Ecotypes,", {"font": HEAD_FONT, "size": 52, "bold": True, "color": WHITE})],
    [("One Genome", {"font": HEAD_FONT, "size": 52, "bold": True, "color": ICE})],
], space_after=2)
txt(s, 0.92, 4.0, 11.5, 0.8,
    [[("How presence-absence variation and DNA methylation may drive ",
       {"size": 21, "color": ICE}),
      ("lean vs. siscowet", {"size": 21, "bold": True, "color": WHITE}),
      (" divergence in lake trout", {"size": 21, "color": ICE})]])
txt(s, 0.92, 6.95, 8, 0.4, [[("15-minute research talk  ·  PacBio HiFi comparative genomics",
                              {"size": 13, "color": WHITE})]])

# ============================================================ 2. THE QUESTION
s = slide(WHITE)
kicker(s, 0.9, 0.6, "The biological question")
txt(s, 0.9, 1.0, 11.5, 1.0,
    [[("What makes a lean a lean and a siscowet a siscowet?",
       {"font": HEAD_FONT, "size": 33, "bold": True, "color": INK})]])

# two ecotype cards
def ecocard(x, color, name, depth, traits):
    rect(s, x, 2.15, 5.4, 3.9, LIGHTBG, shadow=True, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    rect(s, x, 2.15, 5.4, 0.85, color, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    rect(s, x, 2.55, 5.4, 0.45, color)  # square off bottom of header
    txt(s, x+0.4, 2.32, 4.6, 0.5, [[(name, {"font": HEAD_FONT, "size": 24, "bold": True, "color": WHITE})]])
    txt(s, x+0.4, 3.15, 4.6, 0.4, [[(depth, {"size": 15, "bold": True, "italic": True, "color": color})]])
    bullets(s, x+0.4, 3.7, 4.7, 2.2, traits, size=15, gap=9, bullet_color=color)

ecocard(0.9, TEAL, "Lean", "Shallow · nearshore",
        ["Streamlined, low lipid", "Forages shallow water", "Reference-like body plan"])
ecocard(7.05, MID, "Siscowet", "Deep · offshore",
        ["High lipid / fat stores", "Adapted to depth & pressure", "Distinct morphology"])
txt(s, 0.9, 6.35, 11.5, 0.7,
    [[("Same species, same reference genome — so where is the difference written? ",
       {"size": 17, "color": INK}),
      ("We look at two layers.", {"size": 17, "bold": True, "color": CORAL})]])
page_no(s, num())

# ============================================================ 3. TWO-LAYER FRAMEWORK
s = slide(DEEP)
kicker(s, 0.9, 0.6, "The framework", ICE)
txt(s, 0.9, 1.0, 11.5, 0.9, [[("Two layers of genomic difference",
                               {"font": HEAD_FONT, "size": 33, "bold": True, "color": WHITE})]])

rect(s, 0.9, 2.2, 5.5, 3.9, MID, shadow=True, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
txt(s, 1.3, 2.55, 4.8, 0.5, [[("LAYER 1 · PAV", {"size": 15, "bold": True, "color": CORAL, "spacing": 1.5})]])
txt(s, 1.3, 3.0, 4.8, 0.6, [[("The hardware", {"font": HEAD_FONT, "size": 26, "bold": True, "color": WHITE})]])
txt(s, 1.3, 3.75, 4.8, 2.1,
    [[("Presence-absence variation — segments of DNA ", {"size": 16, "color": ICE}),
      ("gained or lost", {"size": 16, "bold": True, "color": WHITE}),
      (" between ecotypes.", {"size": 16, "color": ICE})],
     [("", {})],
     [("Changes ", {"size": 16, "color": ICE}),
      ("what genes are present.", {"size": 16, "bold": True, "color": WHITE})]])

rect(s, 6.9, 2.2, 5.5, 3.9, TEAL, shadow=True, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
txt(s, 7.3, 2.55, 4.8, 0.5, [[("LAYER 2 · 5mC", {"size": 15, "bold": True, "color": RGBColor(0x9B,0xE8,0xDC), "spacing": 1.5})]])
txt(s, 7.3, 3.0, 4.8, 0.6, [[("The software", {"font": HEAD_FONT, "size": 26, "bold": True, "color": WHITE})]])
txt(s, 7.3, 3.75, 4.8, 2.1,
    [[("DNA methylation — reversible marks that ", {"size": 16, "color": ICE}),
      ("switch genes up or down", {"size": 16, "bold": True, "color": WHITE}),
      (".", {"size": 16, "color": ICE})],
     [("", {})],
     [("Changes ", {"size": 16, "color": ICE}),
      ("how genes are used.", {"size": 16, "bold": True, "color": WHITE})]])

txt(s, 0.9, 6.45, 11.5, 0.6,
    [[("Same PacBio HiFi reads resolve both layers at once.",
       {"size": 17, "italic": True, "bold": True, "color": ICE})]], align=PP_ALIGN.CENTER)
page_no(s, num())

# ============================================================ 4. STUDY DESIGN
s = slide(WHITE)
kicker(s, 0.9, 0.6, "Study design & data")
txt(s, 0.9, 1.0, 11.5, 0.9, [[("Eight fish, one platform, two readouts",
                               {"font": HEAD_FONT, "size": 33, "bold": True, "color": INK})]])

# stat tiles
def tile(x, big, small, color):
    rect(s, x, 2.2, 2.75, 1.7, LIGHTBG, shadow=True, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    rect(s, x, 2.2, 0.12, 1.7, color)
    txt(s, x+0.32, 2.42, 2.3, 0.9, [[(big, {"font": HEAD_FONT, "size": 40, "bold": True, "color": color})]])
    txt(s, x+0.34, 3.32, 2.3, 0.5, [[(small, {"size": 13, "color": MUTE})]])

tile(0.9, "8", "PacBio HiFi individuals", MID)
tile(3.95, "4 + 4", "lean + siscowet", TEAL)
tile(7.0, "1", "reference: SaNama_1.0", SEAFOAM)
tile(10.05, "2", "layers from one dataset", CORAL)

rect(s, 0.9, 4.25, 11.5, 1.75, LIGHTBG, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
txt(s, 1.25, 4.5, 11, 0.45, [[("Why HiFi matters here",
                               {"font": HEAD_FONT, "size": 18, "bold": True, "color": INK})]])
bullets(s, 1.25, 5.0, 11, 1.0,
        [[("Long, accurate reads call ", {"size": 15}), ("structural variants (PAV)", {"size": 15, "bold": True}),
          (" and ", {"size": 15}), ("5mC methylation", {"size": 15, "bold": True}),
          (" from the very same molecules — no separate assays, no batch confound.", {"size": 15})]],
        gap=6)
txt(s, 0.9, 6.35, 11.5, 0.6,
    [[("Samples — Lean: bc2041/2068/2069/2070   ·   Siscowet: bc2071/2072/2073/2096   ·   Supporting layer: liver RNA-seq (202 DETs)",
       {"size": 12.5, "color": MUTE})]])
page_no(s, num())

# ============================================================ 5. PAV — what it is
s = slide(WHITE)
kicker(s, 0.9, 0.6, "Layer 1 · Presence-Absence Variation")
txt(s, 0.9, 1.0, 7.5, 1.6, [[("Whole segments of DNA present in one ecotype, absent in the other",
                              {"font": HEAD_FONT, "size": 30, "bold": True, "color": INK})]])
bullets(s, 0.9, 3.0, 7.2, 3.0,
        [[("Detected two ways from HiFi alignments:", {"size": 17, "bold": True})],
         [("Deletions", {"size": 16, "bold": True, "color": MID}), (" — coverage drops to zero over a region", {"size": 16})],
         [("Insertions", {"size": 16, "bold": True, "color": CORAL}), (" — novel sequence flagged in read CIGAR strings", {"size": 16})],
         [("Classified as lean-specific, siscowet-specific, or shared", {"size": 16})]],
        gap=12)
# side visual: two strands, one with a gap
vx = 8.7
rect(s, vx, 2.4, 3.7, 3.6, LIGHTBG, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
txt(s, vx, 2.55, 3.7, 0.4, [[("Lean", {"size": 14, "bold": True, "color": TEAL})]], align=PP_ALIGN.CENTER)
for i,(seg,col) in enumerate([(0,MID),(1,MID),(2,CORAL),(3,MID)]):
    rect(s, vx+0.35+i*0.78, 3.0, 0.68, 0.45, col)
txt(s, vx, 4.05, 3.7, 0.4, [[("Siscowet", {"size": 14, "bold": True, "color": MID})]], align=PP_ALIGN.CENTER)
for i,(seg,col,present) in enumerate([(0,MID,1),(1,RGBColor(0xD7,0xDE,0xE2),0),(2,CORAL,1),(3,MID,1)]):
    rect(s, vx+0.35+i*0.78, 4.45, 0.68, 0.45, col)
txt(s, vx+0.35+0.78, 4.95, 0.9, 0.5, [[("deletion", {"size": 11, "italic": True, "color": MUTE})]], align=PP_ALIGN.CENTER)
txt(s, vx, 5.5, 3.7, 0.4, [[("same locus → different content", {"size": 12, "italic": True, "color": MUTE})]], align=PP_ALIGN.CENTER)
page_no(s, num())

# ============================================================ 6. PAV results (chart)
s = slide(WHITE)
kicker(s, 0.9, 0.6, "Layer 1 · Results")
txt(s, 0.9, 1.0, 11.5, 0.9, [[("Siscowet carries more ecotype-specific structure",
                               {"font": HEAD_FONT, "size": 31, "bold": True, "color": INK})]])

cd = CategoryChartData()
cd.categories = ["Lean-specific", "Siscowet-specific", "Shared"]
cd.add_series("Insertions", (770891, 1086799, 0))
cd.add_series("Deletions", (225337, 245906, 0))
cd.add_series("Shared variants", (0, 0, 878372))
gframe = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_STACKED,
                            Inches(0.9), Inches(2.05), Inches(7.4), Inches(4.4), cd)
chart = gframe.chart
chart.has_title = False
chart.has_legend = True
chart.legend.position = XL_LEGEND_POSITION.BOTTOM
chart.legend.include_in_layout = False
chart.legend.font.size = Pt(12)
plot = chart.plots[0]
plot.gap_width = 80
cols = [CORAL, MID, TEAL]
for ser, c in zip(chart.series, cols):
    ser.format.fill.solid()
    ser.format.fill.fore_color.rgb = c
cat = chart.category_axis
cat.tick_labels.font.size = Pt(12)
cat.tick_labels.font.bold = True
val = chart.value_axis
val.tick_labels.font.size = Pt(10)
val.has_major_gridlines = True

# big number callouts on right
def callout(y, num_s, label, color):
    txt(s, 8.7, y, 3.9, 0.7, [[(num_s, {"font": HEAD_FONT, "size": 33, "bold": True, "color": color})]])
    txt(s, 8.72, y+0.62, 3.9, 0.4, [[(label, {"size": 13, "color": MUTE})]])
callout(2.15, "996,228", "lean-specific variants", CORAL)
callout(3.45, "1,332,705", "siscowet-specific variants", MID)
callout(4.75, "878,372", "shared between ecotypes", TEAL)
txt(s, 8.7, 5.85, 3.9, 0.8, [[("≈ 34% more ecotype-specific variation in siscowet",
                               {"size": 14, "italic": True, "bold": True, "color": INK})]])
page_no(s, num())

# ============================================================ 7. PAV mechanism
s = slide(LIGHTBG)
kicker(s, 0.9, 0.6, "Layer 1 · From variants to phenotype")
txt(s, 0.9, 1.0, 11.5, 0.9, [[("How presence-absence becomes ecotype",
                               {"font": HEAD_FONT, "size": 32, "bold": True, "color": INK})]])
def mech(x, title, body, color, icon):
    rect(s, x, 2.2, 3.65, 3.4, WHITE, shadow=True, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    rect(s, x+0.35, 2.55, 0.7, 0.7, color, shape=MSO_SHAPE.OVAL)
    txt(s, x+0.35, 2.62, 0.7, 0.55, [[(icon, {"size": 24, "bold": True, "color": WHITE})]], align=PP_ALIGN.CENTER)
    txt(s, x+0.35, 3.45, 3.0, 0.6, [[(title, {"font": HEAD_FONT, "size": 19, "bold": True, "color": INK})]])
    txt(s, x+0.35, 4.15, 3.0, 1.3, [[(body, {"size": 14.5, "color": INK})]])
mech(0.9, "Add / remove genes", "A PAV overlapping a coding region adds or deletes a gene outright — a present-vs-absent function.", MID, "±")
mech(4.85, "Shift gene dosage", "Copy-number change tunes expression up or down without touching the coding sequence.", TEAL, "×")
mech(8.8, "Rewire regulation", "Indels in promoters / enhancers reposition regulatory elements, altering control.", CORAL, "~")
txt(s, 0.9, 5.95, 11.5, 1.0,
    [[("Candidate targets for ecotype divergence: ", {"size": 16, "color": INK}),
      ("lipid metabolism, depth & pressure tolerance, sensory adaptation", {"size": 16, "bold": True, "color": MID}),
      (" — flag PAV-hit genes here and test against expression.", {"size": 16, "color": INK})]])
page_no(s, num())

# ============================================================ 8. Methylation intro
s = slide(WHITE)
kicker(s, 0.9, 0.6, "Layer 2 · DNA Methylation", SEAFOAM)
txt(s, 0.9, 1.0, 7.6, 1.6, [[("Same sequence, different settings",
                              {"font": HEAD_FONT, "size": 31, "bold": True, "color": INK})]])
bullets(s, 0.9, 2.9, 7.3, 3.2,
        [[("5-methylcytosine (5mC) at CpG sites", {"size": 16, "bold": True, "color": SEAFOAM}),
          (" — read directly off HiFi data", {"size": 16})],
         [("Methylation typically ", {"size": 16}), ("silences", {"size": 16, "bold": True}),
          (" genes and represses transposable elements", {"size": 16})],
         [("Reversible & potentially heritable", {"size": 16, "bold": True}),
          (" — a route to plastic phenotype", {"size": 16})],
         [("We test each CpG for ecotype differences, then group into regions (DMRs)", {"size": 16})]],
        gap=13, bullet_color=SEAFOAM)
# side: methyl tag concept
vx = 8.8
rect(s, vx, 2.5, 3.6, 3.3, LIGHTBG, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
rect(s, vx+0.4, 4.4, 2.8, 0.5, RGBColor(0xC9,0xD2,0xD7))
for i in range(5):
    on = i in (1,3)
    c = SEAFOAM if on else RGBColor(0x9A,0xA8,0xAF)
    cx = vx+0.55+i*0.55
    rect(s, cx, 4.45, 0.16, 0.4, RGBColor(0x6B,0x7A,0x82))  # CpG tick
    if on:
        rect(s, cx-0.12, 3.95, 0.4, 0.4, c, shape=MSO_SHAPE.OVAL)
        txt(s, cx-0.12, 4.0, 0.4, 0.32, [[("m", {"size": 14, "bold": True, "color": WHITE})]], align=PP_ALIGN.CENTER)
txt(s, vx, 5.15, 3.6, 0.5, [[("methyl marks tune the same DNA", {"size": 12.5, "italic": True, "color": MUTE})]], align=PP_ALIGN.CENTER)
page_no(s, num())

# ============================================================ 9. Methylation results
s = slide(WHITE)
kicker(s, 0.9, 0.6, "Layer 2 · Results", SEAFOAM)
txt(s, 0.9, 1.0, 11.5, 0.9, [[("From half a million sites to 302 regions",
                               {"font": HEAD_FONT, "size": 31, "bold": True, "color": INK})]])
# funnel of three stat tiles
def funnel(x, w, big, label, color):
    rect(s, x, 2.25, w, 1.55, LIGHTBG, shadow=True, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    rect(s, x, 2.25, w, 0.14, color)
    txt(s, x, 2.5, w, 0.8, [[(big, {"font": HEAD_FONT, "size": 34, "bold": True, "color": color})]], align=PP_ALIGN.CENTER)
    txt(s, x, 3.32, w, 0.4, [[(label, {"size": 13, "color": MUTE})]], align=PP_ALIGN.CENTER)
funnel(0.9, 3.4, "540,040", "CpG sites tested", TEAL)
txt(s, 4.35, 2.85, 0.55, 0.5, [[("→", {"size": 26, "bold": True, "color": MUTE})]], align=PP_ALIGN.CENTER)
funnel(4.95, 3.4, "4,440", "significant DMCs (p < 0.05)", MID)
txt(s, 8.4, 2.85, 0.55, 0.5, [[("→", {"size": 26, "bold": True, "color": MUTE})]], align=PP_ALIGN.CENTER)
funnel(9.0, 3.4, "302", "differentially methylated regions", SEAFOAM)

# directional split chart
cd = CategoryChartData()
cd.categories = ["Hypermethylated", "Hypomethylated"]
cd.add_series("DMRs in siscowet", (20, 282))
gframe = s.shapes.add_chart(XL_CHART_TYPE.BAR_CLUSTERED,
                            Inches(0.9), Inches(4.15), Inches(6.3), Inches(2.7), cd)
ch = gframe.chart
ch.has_title = False
ch.has_legend = False
pl = ch.plots[0]
pl.has_data_labels = True
pl.data_labels.font.size = Pt(13)
pl.data_labels.font.bold = True
pl.data_labels.position = XL_LABEL_POSITION.OUTSIDE_END
ser = ch.series[0]
# color points individually
from pptx.oxml.ns import qn
pts = [(SEAFOAM), (SEAFOAM)]
ser.format.fill.solid()
ser.format.fill.fore_color.rgb = SEAFOAM
ch.category_axis.tick_labels.font.size = Pt(12)
ch.category_axis.tick_labels.font.bold = True
ch.value_axis.tick_labels.font.size = Pt(9)

txt(s, 7.5, 4.4, 4.9, 2.4,
    [[("The asymmetry is the story.", {"size": 18, "bold": True, "font": HEAD_FONT, "color": INK})],
     [("", {"size": 6})],
     [("282 of 302 DMRs are ", {"size": 16, "color": INK}),
      ("hypomethylated", {"size": 16, "bold": True, "color": SEAFOAM}),
      (" in siscowet — a genome-wide tilt toward de-repression.", {"size": 16, "color": INK})]])
page_no(s, num())

# ============================================================ 10. Methylation mechanism
s = slide(DEEP)
kicker(s, 0.9, 0.6, "Layer 2 · Interpreting the tilt", ICE)
txt(s, 0.9, 1.0, 11.5, 0.9, [[("Why a hypomethylated siscowet genome matters",
                               {"font": HEAD_FONT, "size": 30, "bold": True, "color": WHITE})]])
def dcard(x, title, body):
    rect(s, x, 2.25, 3.65, 3.5, MID, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    txt(s, x+0.35, 2.55, 3.0, 0.9, [[(title, {"font": HEAD_FONT, "size": 19, "bold": True, "color": WHITE})]])
    txt(s, x+0.35, 3.55, 3.0, 2.0, [[(body, {"size": 14.5, "color": ICE})]])
dcard(0.9, "Genes switch on", "Lower methylation near promoters can de-repress genes — potentially the lipid & depth programs that define siscowet.")
dcard(4.85, "Mobile elements wake", "Genome-wide hypomethylation can re-activate transposable elements, a documented engine of rapid divergence.")
dcard(8.8, "A heritable dial", "Methylation state can persist across generations — divergence without waiting for new mutations.")
txt(s, 0.9, 6.1, 11.5, 0.7,
    [[("Two independent signals, same direction: siscowet shows ", {"size": 16, "color": ICE}),
      ("more structural variation and broad hypomethylation.", {"size": 16, "bold": True, "color": WHITE})]],
    align=PP_ALIGN.CENTER)
page_no(s, num())

# ============================================================ 11. Integration
s = slide(WHITE)
kicker(s, 0.9, 0.6, "Integration · The key move")
txt(s, 0.9, 1.0, 11.5, 0.9, [[("Do the layers converge on the same genes?",
                               {"font": HEAD_FONT, "size": 31, "bold": True, "color": INK})]])
# three overlapping circles (Venn) — sits in lower 2/3, labels inside lobes
r = 1.45
c1c = (4.5, 4.85)              # bottom-left  (Methylation)
c2c = (5.95, 4.85)             # bottom-right (Expression)
c3c = (5.225, 3.6)             # top          (Structure/PAV)
c1 = rect(s, c1c[0]-r, c1c[1]-r, 2*r, 2*r, CORAL, shape=MSO_SHAPE.OVAL)
c2 = rect(s, c2c[0]-r, c2c[1]-r, 2*r, 2*r, SEAFOAM, shape=MSO_SHAPE.OVAL)
c3 = rect(s, c3c[0]-r, c3c[1]-r, 2*r, 2*r, MID, shape=MSO_SHAPE.OVAL)
for sp in (c1, c2, c3):
    sppr = sp.fill._xPr.find(qn('a:solidFill'))
    clr = sppr.find(qn('a:srgbClr'))
    a = clr.makeelement(qn('a:alpha'), {'val': '58000'})
    clr.append(a)
# labels inside the non-overlapping part of each lobe
txt(s, c3c[0]-1.1, c3c[1]-1.15, 2.2, 0.4, [[("Structure (PAV)", {"size": 13, "bold": True, "color": WHITE})]], align=PP_ALIGN.CENTER)
txt(s, c1c[0]-1.25, c1c[1]+0.45, 1.7, 0.4, [[("Methylation", {"size": 13, "bold": True, "color": WHITE})]], align=PP_ALIGN.CENTER)
txt(s, c2c[0]-0.45, c2c[1]+0.45, 1.7, 0.4, [[("Expression", {"size": 13, "bold": True, "color": WHITE})]], align=PP_ALIGN.CENTER)
txt(s, c3c[0]-1.0, 4.5, 2.0, 0.6, [[("causal", {"size": 12, "bold": True, "color": INK})], [("candidates", {"size": 12, "bold": True, "color": INK})]], align=PP_ALIGN.CENTER)

bullets(s, 8.2, 2.3, 4.4, 4.0,
        [[("Are DMRs enriched near PAV breakpoints?", {"size": 15.5})],
         [("Do PAV-hit genes overlap the 202 liver DETs?", {"size": 15.5})],
         [("Three-way agreement", {"size": 15.5, "bold": True, "color": MID}),
          (" turns correlation into mechanism", {"size": 15.5})],
         [("Next: ", {"size": 15.5, "bold": True}),
          ("ecotype-specific hifiasm assemblies remove reference bias on PAV calls", {"size": 15.5})]],
        gap=14, bullet_color=MID)
page_no(s, num())

# ============================================================ 12. Takeaways
s = slide(LIGHTBG)
kicker(s, 0.9, 0.6, "Take-home")
txt(s, 0.9, 1.0, 11.5, 0.9, [[("Three things to remember",
                               {"font": HEAD_FONT, "size": 33, "bold": True, "color": INK})]])
def take(y, n, head, body, color):
    rect(s, 0.9, y, 11.5, 1.45, WHITE, shadow=True, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    rect(s, 0.9, y, 0.14, 1.45, color)
    rect(s, 1.25, y+0.32, 0.8, 0.8, color, shape=MSO_SHAPE.OVAL)
    txt(s, 1.25, y+0.4, 0.8, 0.6, [[(n, {"font": HEAD_FONT, "size": 26, "bold": True, "color": WHITE})]], align=PP_ALIGN.CENTER)
    txt(s, 2.4, y+0.22, 9.6, 0.5, [[(head, {"font": HEAD_FONT, "size": 19, "bold": True, "color": INK})]])
    txt(s, 2.4, y+0.78, 9.6, 0.55, [[(body, {"size": 14.5, "color": MUTE})]])
take(2.15, "1", "Ecotype divergence is written at two levels",
     "Structural (what genes exist) and epigenetic (how they're used) — not either/or.", MID)
take(3.78, "2", "Both signals point the same way",
     "Siscowet shows ~34% more ecotype-specific PAV and 282/302 DMRs hypomethylated.", CORAL)
take(5.41, "3", "Integration is the next experiment",
     "Overlap PAV, methylation, and expression to find true causal candidates.", SEAFOAM)
page_no(s, num())

# ============================================================ 13. Closing
s = slide(DEEP)
rect(s, 0, 5.4, 13.333, 2.1, MID)
rect(s, 0, 6.4, 13.333, 1.1, TEAL)
rect(s, 0, 7.05, 13.333, 0.45, SEAFOAM)
txt(s, 0.9, 1.7, 11.5, 1.6, [
    [("Two layers, one divergence", {"font": HEAD_FONT, "size": 40, "bold": True, "color": WHITE})],
])
txt(s, 0.92, 3.1, 11, 1.2,
    [[("Presence-absence variation and DNA methylation each track the lean–siscowet split — ",
       {"size": 19, "color": ICE})],
     [("and where they converge, we find the genes that make an ecotype.",
       {"size": 19, "bold": True, "color": WHITE})]], space_after=4)
txt(s, 0.92, 5.7, 11.5, 0.5,
    [[("Explore the data — interactive genome browser:", {"size": 14, "color": ICE})]])
txt(s, 0.92, 6.6, 11.5, 0.5,
    [[("sr320.github.io/project-lake-trout/genome-browser", {"size": 16, "bold": True, "color": WHITE})]])
txt(s, 0.92, 7.12, 8, 0.35, [[("Roberts Lab · University of Washington · School of Aquatic & Fishery Sciences",
                               {"size": 11.5, "color": ICE})]])

out = "/Users/sr320/GitHub/project-lake-trout/output/deck/lake-trout-ecotypes-PAV-methylation.pptx"
prs.save(out)
print("Saved", out, "slides:", len(prs.slides._sldIdLst))
