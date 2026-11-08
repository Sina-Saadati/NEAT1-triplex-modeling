import os
from pymol import cmd, util
D = os.path.abspath("data"); OUT = os.path.abspath("figures"); os.makedirs(OUT, exist_ok=True)
for n, c in [("c_u", [0.902, 0.624, 0.0]), ("c_a", [0.835, 0.369, 0.0]), ("c_g", [0.78, 0.78, 0.78]),
             ("c_p", [0.80, 0.475, 0.655]), ("c_q", [0.94, 0.89, 0.26])]:
    cmd.set_color(n, c)
cmd.bg_color("white")
for k, v in [("antialias", 2), ("ray_shadow", 0), ("depth_cue", 0), ("ray_trace_mode", 1), ("cartoon_ring_mode", 3),
             ("cartoon_ring_finder", 1), ("cartoon_tube_radius", 0.6), ("cartoon_nucleic_acid_mode", 4),
             ("ray_opaque_background", 1), ("specular", 0.3), ("orthoscopic", 1), ("stick_radius", 0.22)]:
    cmd.set(k, v)

def base(name):
    cmd.hide("everything", name); cmd.show("cartoon", name); cmd.color("c_g", name)
    cmd.color("c_u", name + " and resi 14-23"); cmd.color("c_a", name + " and resi 82-101")
def fit(sel, buf):
    cmd.orient(sel); cmd.zoom(sel, buffer=buf, complete=1); cmd.clip("slab", 600)
def shot(png, w=3200, h=2200):
    cmd.png(os.path.join(OUT, png), width=w, height=h, dpi=300, ray=1)

# Four predicted models
for name, fn, png in [("rf_nomsa", "NEAT1_ENE_A_RhoFold.pdb", "3D_RhoFold_noMSA.png"),
                      ("rf_msa", "NEAT1_ENE_A_RhoFold_withMSA.pdb", "3D_RhoFold_withMSA.png"),
                      ("trr", "model1.pdb", "3D_trRosettaRNA.png"),
                      ("rnac", "NEAT1_ENE_A_RNAComposer.pdb", "3D_RNAComposer.png")]:
    p = os.path.join(D, fn)
    if not os.path.exists(p): print("MISSING:", p); continue
    cmd.delete("all"); cmd.load(p, name); base(name); fit(name, 10); shot(png)

# MD frame (100 ps) + pocket 2 + docked quercetin (pose 1)
cmd.delete("all")
cmd.load(os.path.join(D, "md_final_RNA_only.pdb"), "rna"); cmd.load(os.path.join(D, "pocket2_atm.pdb"), "pocket")
cmd.load(os.path.join(D, "quercetin_docked.pdbqt"), "que_all"); cmd.create("que", "que_all", 1, 1); cmd.delete("que_all")
base("rna")
res = "rna and resi 7+10+86+88+89+90"
cmd.show("sticks", res); cmd.color("c_p", res); util.cnc(res)
cmd.hide("everything", "pocket"); cmd.show("surface", "pocket"); cmd.color("c_p", "pocket"); cmd.set("transparency", 0.5, "pocket")
cmd.show("sticks", "que"); cmd.color("c_q", "que"); util.cnc("que")
fit("rna", 10); shot("3D_MD_pocket_overview.png")
cmd.select("near_rna", "byres (rna within 16 of (pocket or que))"); cmd.hide("everything", "rna and not near_rna")
cmd.orient("near_rna or pocket or que"); cmd.zoom("near_rna or pocket or que", buffer=3, complete=1); cmd.clip("slab", 600)
shot("3D_MD_pocket_closeup.png", 2600, 2600)
print("Done ->", OUT)