#!/usr/bin/env python3
"""Workshop talk deck (27 Nov 2026), built on the user's previous deck (14h00_PLT-2551_Loewe.pptx) so the master, layouts (orange header), title slide
and logos are exactly his: its title slide is kept and edited, all other slides are dropped (their media are not written), and 10 content slides are
added on the 'Quote' layout. Figures: assets/clean/*.png (crop_panels.py + clean_panels.py; talk_figs.py for the new graphics); movies: paper SI movies
S1-S4 with poster frames. Speaker notes = the script in PROPOSAL.md.
Usage: python build_deck.py <old deck .pptx> -> Loewe_workshop_20261127.pptx (python-pptx)"""
import os, sys, copy
from PIL import Image
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
H = os.path.dirname(os.path.abspath(__file__)); C = os.path.join(H, "assets", "clean"); A = os.path.join(H, "assets")
ORANGE, INK, GREY, GREEN = RGBColor(0xF6, 0x70, 0x31), RGBColor(0x1A, 0x1A, 0x1A), RGBColor(0x59, 0x59, 0x59), RGBColor(0x1A, 0x7A, 0x3A)
CITE = "B. Loewe & T. N. Shendruk, in preparation"

prs = Presentation(sys.argv[1])
ids = prs.slides._sldIdLst
for sid in list(ids)[1:]:                                                  # keep only the title slide; dropped slides' media are not saved
    prs.part.drop_rel(sid.rId); ids.remove(sid)
QUOTE = next(l for l in prs.slide_layouts if l.name == "Quote")


def fit(slide, path, x, y, w, h, align="center"):
    """picture scaled to fit the box (x, y, w, h in inches), centred (or left/top aligned)"""
    iw, ih = Image.open(path).size; s = min(w / iw, h / ih); pw, ph = iw * s, ih * s
    px = x + (w - pw) / 2 if align == "center" else x; py = y + (h - ph) / 2 if align == "center" else y
    return slide.shapes.add_picture(path, Inches(px), Inches(py), Inches(pw), Inches(ph))


def movie(slide, name, x, y, w, h):
    poster = os.path.join(A, name + "_poster.png"); iw, ih = Image.open(poster).size; s = min(w / iw, h / ih); pw, ph = iw * s, ih * s
    return slide.shapes.add_movie(os.path.join(A, name + ".mp4"), Inches(x + (w - pw) / 2), Inches(y + (h - ph) / 2), Inches(pw), Inches(ph),
                                  poster_frame_image=poster, mime_type="video/mp4")


def text(slide, x, y, w, h, paras, size=18, color=INK, bold=False, anchor=MSO_ANCHOR.TOP, font="Arial", align=PP_ALIGN.LEFT, space=6):
    """paras: list of str or list of (str, dict) runs; a str starting with '• ' becomes a bulleted-looking paragraph with hanging indent"""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)); tf = tb.text_frame; tf.word_wrap = True; tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = Inches(0.05); tf.margin_top = tf.margin_bottom = Inches(0.03)
    for i, para in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph(); p.alignment = align; p.space_after = Pt(space)
        runs = [(para, {})] if isinstance(para, str) else para
        for s, o in runs:
            r = p.add_run(); r.text = s; f = r.font; f.name = o.get("font", font); f.size = Pt(o.get("size", size)); f.bold = o.get("bold", bold)
            f.italic = o.get("italic", False); f.color.rgb = o.get("color", color)
    return tb


def new(title, notes, cite=True):
    s = prs.slides.add_slide(QUOTE)
    for ph in list(s.placeholders): ph._element.getparent().remove(ph._element)   # layout body / slide-number placeholders: not used
    text(s, 0.2, 0.06, 12.9, 1.0, [title], size=26, color=RGBColor(0xFF, 0xFF, 0xFF), font="Nunito Sans", anchor=MSO_ANCHOR.MIDDLE, space=0)
    if cite: text(s, 8.3, 7.08, 4.95, 0.34, [CITE], size=12, color=GREY, align=PP_ALIGN.RIGHT, space=0)
    s.notes_slide.notes_text_frame.text = notes
    return s


# ---------- 1. title slide (kept from the old deck, edited) ----------
t = prs.slides[0]
def set_text(shape, lines):
    tf = shape.text_frame; ps = tf.paragraphs
    for p in ps[len(lines):]: p._p.getparent().remove(p._p)
    for p, line in zip(tf.paragraphs, lines):
        rs = p.runs
        for r in rs[1:]: r._r.getparent().remove(r._r)
        rs[0].text = line
    while len(tf.paragraphs) < len(lines):
        np_ = copy.deepcopy(tf.paragraphs[-1]._p); tf._txBody.append(np_)
        tf.paragraphs[-1].runs[0].text = lines[len(tf.paragraphs) - 1]
for sh in list(t.shapes):
    if sh.name == "Title 1": set_text(sh, ["Geometry, not topology, places the charge of active nematic cavities"])
    elif sh.name == "TextBox 4":                                       # presenter line keeps the old style (underlined); co-author plain, no old line breaks
        tf = sh.text_frame
        for p_ in tf.paragraphs[1:]: p_._p.getparent().remove(p_._p)
        for br in tf.paragraphs[0]._p.findall("{http://schemas.openxmlformats.org/drawingml/2006/main}br"): tf.paragraphs[0]._p.remove(br)
        rs = tf.paragraphs[0].runs
        for r in rs[1:]: r._r.getparent().remove(r._r)
        rs[0].text = "Benjamin Loewe (PUC)"
        np_ = copy.deepcopy(tf.paragraphs[0]._p); tf._txBody.append(np_); r2 = tf.paragraphs[1].runs[0]; r2.text = "Tyler N. Shendruk (Edinburgh)"; r2.font.underline = False
    elif sh.name == "Picture 6": sh._element.getparent().remove(sh._element)
fit(t, os.path.join(C, "f1d.png"), 7.15, 0.6, 5.7, 6.3)
t.notes_slide.notes_text_frame.text = ("Thank you. I'll talk about what happens to topological charge when you confine an active nematic in a cavity whose walls force a "
                                      "winding. This is work with Tyler Shendruk in Edinburgh.")

# ---------- 1b. intro: active nematics and cavities (user, 23 Nov: softer intro for an active-matter audience) ----------
s = new("Active nematics make their own defects; anchored walls impose them",
        "Let me start with the two ingredients. An active nematic is made of elongated units, here extensile: each pushes fluid out along its axis. This active stress makes the "
        "material flow on its own, and in the bulk the flow becomes chaotic, what we call active turbulence. The orientation field is full of topological defects of charge plus "
        "and minus one half, created and annihilated in pairs. The plus one-half has a comet shape and swims head first; the minus one-half is symmetric and does not swim.\n"
        "Now put it in a cavity whose wall fixes the orientation. With homeotropic anchoring, perpendicular to the wall, the director turns once as you go around, so the cavity "
        "must hold a total charge of plus one. A passive nematic splits it into two plus one-half defects that sit still. In the active one the pair keeps moving around the cavity. "
        "Same charge, very different dynamics. The question of this talk is what happens when the wall forces more charge into the cavity.")
text(s, 0.3, 1.3, 6.1, 0.45, [[("Active nematics", {"bold": True, "color": ORANGE})]], size=20, space=0)
fit(s, os.path.join(C, "rods.png"), 0.3, 1.8, 3.05, 0.85)
fit(s, os.path.join(C, "defects.png"), 0.25, 2.7, 3.15, 2.15)
movie(s, "bulk_movie", 3.45, 1.8, 3.0, 2.95)
text(s, 0.3, 4.9, 6.15, 2.0, ["Extensile active stress drives spontaneous flows: active turbulence",
                              "±½ defects are created and annihilated in pairs",
                              "+½ swim head first; −½ do not"], size=16, space=6)
ln = s.shapes.add_connector(1, Inches(6.62), Inches(1.45), Inches(6.62), Inches(6.9)); ln.line.color.rgb = RGBColor(0xBF, 0xBF, 0xBF); ln.line.width = Pt(1.5)
text(s, 6.8, 1.3, 6.3, 0.45, [[("Cavities: the wall fixes the charge", {"bold": True, "color": ORANGE})]], size=20, space=0)
text(s, 6.8, 1.8, 3.0, 0.35, ["Passive"], size=16, align=PP_ALIGN.CENTER, space=0)
text(s, 9.95, 1.8, 3.0, 0.35, ["Active"], size=16, align=PP_ALIGN.CENTER, space=0)
fit(s, os.path.join(A, "k1_passive.png"), 6.8, 2.2, 3.0, 2.75)
movie(s, "k1_active_movie", 9.95, 2.2, 3.0, 2.75)
text(s, 6.8, 5.05, 6.3, 1.9, ["Homeotropic wall: the director turns once (k = 1), so the cavity holds charge +1",
                              "Passive: two +½ at rest; active: the pair keeps moving",
                              [("Same charge, different dynamics. What if the wall imposes more?", {"bold": True})]], size=16, space=6)

# ---------- 2. the question ----------
s = new("Winding k: the cavity must hold charge k, but topology does not say where",
        "Now let the anchoring turn the director k times around the wall. Topology tells you that the cavity contains a total charge k. It does not tell you where that charge sits.\n"
        "For a passive nematic the answer is simple. All of it stays at the wall, as a ring of 2k half-integer defects, and the interior is defect-free.\n"
        "In active nematics, previous work found that the boundary charge is screened: in turbulent disks it stays in a thin "
        "boundary layer, and below the pair-creation threshold the extra defects stay close to the wall.\n"
        "We pushed to large windings, up to k equal to ten, in the regime of continual pair creation. The question is simple: where does the charge go?")
fit(s, os.path.join(C, "anchoring.png"), 0.4, 1.45, 4.0, 4.5)
fit(s, os.path.join(C, "f1a_passive.png"), 4.7, 1.45, 3.3, 4.5)
text(s, 0.4, 6.0, 3.9, 0.9, ["Director turns k times along the wall"], size=16, color=GREY, align=PP_ALIGN.CENTER)
text(s, 4.6, 6.0, 3.5, 0.9, ["Passive nematic: all 2k = 12 +½ at the wall (k = 6)"], size=16, color=GREY, align=PP_ALIGN.CENTER)
text(s, 8.5, 1.6, 4.5, 4.4, [[("Active nematics in disks", {"bold": True})],
                              "Norton et al.: turbulent disks screen the boundary charge in a thin boundary layer",
                              "Klein et al.: windings up to q = 5 below pair creation; defects braid near the wall",
                              [("Here: k up to 10, with continual pair creation", {"bold": True})]], size=18, space=14)
text(s, 8.5, 5.9, 4.5, 0.9, [[("Topology says how much charge, not where.", {"bold": True, "color": ORANGE})]], size=20)

# ---------- 3. setup ----------
s = new("We simulate active nematics in cavities with imposed winding",
        "We use multi-particle collision dynamics, a particle-based solver for nematohydrodynamics, with extensile active dipoles. The relevant scale is the active length, "
        "about twenty cells here. The walls are no-slip and impose strong anchoring that winds k times. We vary the winding, the size from a fifth to three active lengths, "
        "and the shape, from circles to elongated ellipses and ovals. Everything I show is checked against a continuum Q-tensor model with backflow.")
fit(s, os.path.join(C, "f1b.png"), 0.5, 1.45, 5.2, 5.4)
text(s, 6.2, 1.6, 6.7, 5.2, [[("NIAMH-MPCD", {"bold": True}), (": particle-based 2D nematohydrodynamics, extensile active dipoles", {})],
                              "Activity ζ = 0.07: active length ℓa ≈ 22 cells",
                              "Flow aligning: λ = 2 (varied 0.5–6)",
                              "No-slip walls, strong anchoring θ = kφ (winding k)",
                              "Circles: R = 0.2–2.9 ℓa, k = 2–10",
                              "Ellipses and Cassini ovals: aspect 1.1–2.5, three sizes (675 cavities)",
                              [("Check: ", {"bold": True}), ("continuum Q-tensor model with backflow", {})]], size=20, space=12)

# ---------- 4. observation ----------
s = new("An active cavity keeps exactly 2k−2 defects at the wall and sends +1 into the bulk",
        "Here is what an active cavity does, for k equal to six. A passive cavity would have twelve defects at the wall. The active one keeps ten, sharply localised, and these sites "
        "almost never empty. In the bulk, pairs are created and annihilated all the time, but the net charge in the bulk is exactly plus one.\n"
        "You can also see that a single vortex fills the cavity, and the wall defects are all tilted in the same sense, following the vortex.\n"
        "This is not specific to k equal to six. Over windings from two to ten and sizes from about one to two active lengths, the wall holds two k minus two defects and the bulk "
        "plus one. Topology is satisfied, but the boundary does not screen the charge: exactly one unit always moves into the bulk. Too small and the vortex cannot form, too big "
        "and the cavity becomes turbulent.\nSo why two k minus two, and why plus one?")
movie(s, "S1_circle_overcharged", 0.3, 1.3, 5.3, 5.65)
fit(s, os.path.join(C, "f1a_all.png"), 5.8, 1.3, 7.3, 2.9)
fit(s, os.path.join(C, "f1c.png"), 5.8, 4.25, 4.6, 2.75)
text(s, 10.5, 4.35, 2.65, 2.6, [[("k = 2–10", {"bold": True})], [("R = 0.8–2.2 ℓa", {"bold": True})], "Wall occupancy 0.997", "Bulk charge 1.01"], size=16, space=8)

# ---------- 5. geometry ----------
s = new("The anchoring ties each wall defect's tilt to its position",
        "The answer is in the geometry of the anchoring. Relative to the local wall normal, the anchored director turns by (k minus one) times the angle. So it is exactly homeotropic "
        "at two k minus two points around the wall.\n"
        "Now take a wall defect. Its orientation, the direction of its tail, is fixed by where it sits: the tilt from the normal is (two k minus two) times its position, plus a "
        "constant. The simulations follow this lock almost perfectly. So position and tilt are one degree of freedom.\n"
        "That has a consequence. A ring of two k minus two defects can sit on the homeotropic sites, untilted, and it can tilt as a whole simply by sliding all its defects together "
        "along the wall. A ring of two k cannot: with that spacing, the lock makes the tilts wind twice around the ring. There is no way to give all of them the same tilt.\n"
        "In a passive cavity this does not matter, because nothing asks for a common tilt, and the two k ring is perfectly stable.")
fit(s, os.path.join(C, "rings.png"), 0.3, 1.4, 8.9, 5.5)
text(s, 9.5, 1.6, 3.6, 5.2, ["Relative to the wall normal the anchoring turns by (k−1)φ: homeotropic at 2k−2 sites",
                             "Simulations obey the lock: mean resultant length 0.99",
                             [("Passive cavity: nothing asks for a common tilt, so the 2k ring is stable", {"color": GREY})]], size=18, space=16)

# ---------- 6. the flow does it ----------
s = new("A coherent vortex imposes a common tilt, so the 2k ring must eject two",
        "What asks for a common tilt is the flow. The coherent vortex shears the wall layer in one direction everywhere, and flow alignment tilts every wall defect the same way. "
        "For the ring of two k, forcing a common tilt costs roughly twice the elastic energy needed to move two defects into the bulk. So the ring ejects two. What is left is the "
        "largest ring that can tilt as a whole, and the bulk receives plus one.\n"
        "To show that it is the flow and not activity itself, we take a passive cavity and impose a vortex by hand. Once the flow is strong enough to give all wall defects the same "
        "tilt, the ring of twelve ejects exactly two, here in the movie. With weak flow, or with a tumbling nematic that has no preferred angle, nothing happens.\n"
        "And an active continuum model, started from the passive ring, organises itself into exactly the same state.")
movie(s, "S2_passive_swirl_ejection", 0.3, 1.3, 5.3, 5.65)
fit(s, os.path.join(C, "f2c.png"), 5.8, 1.3, 7.3, 4.6)
text(s, 5.9, 5.95, 7.2, 1.1, ["Common tilt on the 2k ring costs ≈ 2× the energy of moving two +½ into the bulk",
                              "Active continuum model, started from the passive ring: same state"], size=16, space=4)

# ---------- 7. tilt law (revised 23 Nov: why twice the Leslie angle) ----------
s = new("The tilt is the flow-aligned angle at the no-slip wall: τ = arccos(1/λ)",
        "Which tilt does the flow select? The key is what a plus one-half defect looks like next to the wall. Its tail points at the wall, where the director is perpendicular to it, "
        "and straight above it, in the bulk, the director lies parallel to the wall. The defect is what turns the director by ninety degrees.\n"
        "Now switch on the shear. The bulk director wants the flow-aligning, Leslie, angle, so everything around the defect rotates by theta L.\n"
        "And here is the geometry: around a plus one-half, the director turns only half as fast as you go around the core. The tail sits where the director points at the core. "
        "If all the directors rotate by theta L, you have to go twice as far around before they point at the core again. So the comet turns by twice the Leslie angle, which is "
        "arccos of one over lambda.\n"
        "The wall lets this happen for free: sliding the defect along the winding anchoring rotates its surroundings at no elastic cost. That's the lock, with the same factor of two.\n"
        "The prediction contains neither the shear rate nor the elastic constants, and the active simulations follow it within three degrees for lambda from two to six.\n\n"
        "[If asked about the angle at the wall vs the bulk] Only the ratio of strain rate to vorticity enters, and the no-slip wall fixes it. In the passive continuum cavity, removing "
        "the shear next to the wall moves the tilt away from arccos(1/lambda) (67 to 96 degrees at lambda 2), while confining all shear to the wall band leaves it unchanged (68).")
fit(s, os.path.join(C, "tilt_factor2.png"), 0.25, 1.35, 8.2, 3.1)
text(s, 0.4, 4.6, 7.8, 2.4, ["Around a +½ the director turns half as fast as you go around the core",
                             "Flow rotates the surroundings by θL, so the comet turns by 2θL = arccos(1/λ)",
                             [("Sliding along the winding anchoring makes this free (the lock, same factor 2)", {"color": GREY})]], size=18, space=10)
fit(s, os.path.join(C, "f3b.png"), 8.6, 1.4, 4.5, 4.5)
text(s, 8.65, 6.0, 4.5, 1.0, [[("No shear rate, no elastic constants: ", {"bold": True}), ("MPCD within 3° for λ = 2–6", {})]], size=16)

# ---------- 8. beyond circles ----------
s = new("Elongating the cavity breaks the vortex, and the wall layer leaks",
        "What if the cavity is not a circle? When we elongate it, the single vortex breaks, at an aspect ratio of about 1.3, for ellipses and ovals and at every size. Without the "
        "vortex, there is no common tilt, and the wall layer starts to lose defects, one at a time. The bulk now holds two, three or more.\n"
        "The released charge is not random. In elongated cavities it follows a wall law: the bulk charge per unit wall length depends only on the winding density, not on the shape "
        "or size, and a straight periodic channel is its flat limit.\n"
        "Stronger flow alignment does not save the layer, which tells us the protection really came from the coherent vortex.")
movie(s, "S3_ellipse_leaking", 0.3, 1.3, 6.3, 5.65)
fit(s, os.path.join(C, "f4a.png"), 6.8, 1.3, 3.05, 4.3)
fit(s, os.path.join(C, "f4b.png"), 10.0, 1.3, 3.1, 4.3)
text(s, 6.85, 5.75, 6.3, 1.25, ["Single vortex breaks at aspect ≈ 1.3 (every size, ellipses and ovals)",
                                "Released charge follows a wall law: bulk charge per wall length set by the winding density"], size=16, space=4)

# ---------- 9. why it leaks ----------
s = new("Untilted wall defects propel themselves out; −½ layers, without jets, do not leak",
        "Why do the defects leave? A plus one-half defect is self-propelled along its head. In the circle the tilt turns this jet along the wall, and the no-slip wall absorbs most "
        "of the inward part. Following single escapes, we see that a defect first untilts over about forty time units. Its head turns into the bulk, and then it runs away, "
        "against an elastic barrier that has not changed.\n"
        "The cleanest test is a mirror control. If we reverse the winding, the same geometry puts minus one-half defects at the wall, two more instead of two fewer. Minus one-half "
        "defects are not self-propelled. In exactly the same elongated cavities, the vortex breaks just as before, but the minus one-half layer does not leak: its escape rate is "
        "more than twenty times lower.\nSo losing the vortex is not enough. The leak needs the wall defects' own jets.")
movie(s, "S4_escape_clean", 0.3, 1.3, 6.6, 4.6)
fit(s, os.path.join(C, "f4e.png"), 7.1, 1.3, 2.95, 3.9)
fit(s, os.path.join(C, "f4g.png"), 10.2, 1.3, 2.9, 3.9)
text(s, 0.4, 6.0, 6.5, 1.0, ["Tilt turns each +½'s jet along the wall; the no-slip wall passes only 30–45 % of its inward part"], size=16, color=GREY)
text(s, 7.15, 5.3, 6.0, 1.7, ["Escapers untilt over ≈ 40 t.u., then leave against an unchanged elastic barrier",
                              [("Mirror (−k): wall −½, no jets: escape rate > 20× lower", {"bold": True})]], size=16, space=6)

# ---------- 10. work in progress: crowding ----------
s = new("Work in progress: crowded walls leak ten times less, and their defects resist the same inward flow",
        "One open question. Densely packed walls leak much less: in aspect-two ellipses, about ten times less per defect. We first thought crowding weakens the inward push. It "
        "doesn't. Measuring the flow at each wall defect, crowded defects get the same inward flow as sparse ones, or more, but they barely follow it: at the same distance from "
        "the wall, their response is about seven times weaker. So crowding strengthens the hold on each defect rather than weakening the drive. Interestingly, the continuum model "
        "does not reproduce this regime: there, crowded layers either freeze or leak freely.\n\n"
        "[If asked] 5 seeds per density, one shape and activity; MPCD flow noise is controlled (equal noise fraction at the compared depth); 'hold' = motion not carried by the "
        "measured flow, elastic in the natural reading, but not yet separated from self-propulsion.")
fit(s, os.path.join(C, "f4b.png"), 0.3, 1.35, 3.6, 4.5)
fit(s, os.path.join(C, "coreflow.png"), 4.1, 1.35, 9.0, 3.9)
text(s, 0.35, 5.95, 3.6, 1.0, ["Vacant fraction falls with wall density ρ"], size=16, color=GREY, align=PP_ALIGN.CENTER)
text(s, 4.2, 5.4, 8.9, 1.6, [[("Same or more inward flow, ≈ 7× weaker response: ", {"bold": True}), ("crowding strengthens the hold on each defect, not a weaker drive", {})],
                             [("Continuum model: crowded layers freeze or leak freely, never this regime", {"color": GREY})]], size=17, space=6)

# ---------- 11. conclusions ----------
s = new("Conclusions",
        "To summarise. Topology fixes how much charge an active cavity holds, but geometry decides where it goes. The wall keeps the largest ring that can tilt as a whole, two k minus "
        "two, and the coherent vortex is what makes that constraint bite, at the flow-aligned tilt. Without the vortex, the wall defects' jets turn inward and the layer leaks, "
        "while a layer of minus one-half defects does not. All of this should be testable in active nematics in wells with patterned anchoring. Thank you.")
text(s, 0.5, 1.5, 7.6, 4.8, ["• Topology fixes how much charge; geometry decides where: the wall keeps the largest ring that can tilt as a whole, 2|k−1|",
                             "• The coherent vortex activates the geometry by imposing a common tilt, τ = arccos(1/λ)",
                             "• Without the vortex (elongated cavities) the jets turn inward and the layer leaks, to a wall law; −½ layers do not",
                             "• Testable: microtubule–kinesin or cell layers in wells with patterned anchoring"], size=20, space=16)
fit(s, os.path.join(C, "f1d.png"), 8.5, 1.4, 4.5, 4.9)
text(s, 0.5, 6.45, 8.0, 0.5, ["Thanks: Tyler N. Shendruk (Edinburgh) · FONDECYT 11250894 · NIAMH-MPCD (Shendruk Lab)"], size=14, color=GREY)

out = os.path.join(H, "Loewe_workshop_20261127.pptx"); prs.save(out); print("wrote", out, os.path.getsize(out) // 1024, "kB")
