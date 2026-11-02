#!/usr/bin/env python3
"""RingLampe – prozedurale Erzeugung (STL, GLB, technische Zeichnung PDF).

Geometrie (Einheiten mm):
  * 2 ringförmige Basen, AD 110 / ID 100, Dicke PLATE_T, mit je 96 viereckigen LED-Öffnungen
  * obere Basis 200 mm über der unteren, um 120° verdreht
  * 96 Röhren, linear gelofted vom unteren zum oberen Segment, Wandstärke 0.2 mm (normal zur Wand)
  * Jeder Horizontalschnitt des Röhrenbündels ist eine Ähnlichkeitsabbildung (Drehung+Skalierung)
    der Basis-Segmente -> Röhren schneiden sich nie, solange die Basis-Segmente disjunkt sind.
"""
import json
import struct
from datetime import date
from pathlib import Path

import numpy as np

OUT = Path(__file__).parent

# ---------------- Parameter ----------------
R_IN, R_OUT = 50.0, 55.0          # Ring ID 100 / AD 110
N_SEG = 96                        # LEDs / Segmente
HEIGHT = 200.0                    # Abstand der Basen
TWIST = np.radians(120.0)         # Verdrehung obere Basis
WALL = 0.2                        # Röhrenwandstärke
PLATE_T = 2.0                     # Dicke der Basisringe
SEG_R_IN, SEG_R_OUT = 50.5, 54.5  # radiale Lage der Segmente (0.5 mm Steg innen/außen)
GAP = 0.6                         # tangentialer Spalt zwischen Segmenten am Innenradius (Basis)
N_SLICES = 80                     # Unterteilung der Röhren in z

PITCH = 2 * np.pi / N_SEG
SEG_W = PITCH - GAP / SEG_R_IN    # Winkelbreite eines Segments


def rot(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s], [s, c]])


def seg_quad(k, twist=0.0):
    """Außenkontur (Footprint) von Segment k, CCW: in-a, out-a, out-b, in-b."""
    th = k * PITCH + twist
    a, b = th - SEG_W / 2, th + SEG_W / 2
    return np.array([
        [SEG_R_IN * np.cos(a), SEG_R_IN * np.sin(a)],
        [SEG_R_OUT * np.cos(a), SEG_R_OUT * np.sin(a)],
        [SEG_R_OUT * np.cos(b), SEG_R_OUT * np.sin(b)],
        [SEG_R_IN * np.cos(b), SEG_R_IN * np.sin(b)],
    ])


def line_intersect(p1, d1, p2, d2):
    A = np.array([d1, -d2]).T
    s = np.linalg.solve(A, p2 - p1)
    return p1 + s[0] * d1


def tube_slices(k):
    """Liefert outer[N+1,4,3], inner[N+1,4,3] für Röhre k."""
    P = seg_quad(k)
    RP = P @ rot(TWIST).T
    ts = np.linspace(0, 1, N_SLICES + 1)
    outer = np.zeros((len(ts), 4, 3))
    inner = np.zeros((len(ts), 4, 3))
    for j, t in enumerate(ts):
        C = (1 - t) * P + t * RP
        z = HEIGHT * t
        cen = C.mean(axis=0)
        lines = []
        for i in range(4):
            c0, c1 = C[i], C[(i + 1) % 4]
            e = np.append(c1 - c0, 0.0)
            # Richtung der Erzeugenden an Kantenmitte
            dmid = ((RP[i] - P[i]) + (RP[(i + 1) % 4] - P[(i + 1) % 4])) / 2
            d = np.append(dmid, HEIGHT)
            n = np.cross(e, d)
            n /= np.linalg.norm(n)
            h = np.hypot(n[0], n[1])
            nh = n[:2] / h
            mid = (c0 + c1) / 2
            if np.dot(cen - mid, nh) < 0:
                nh = -nh
            off = WALL / h                         # horizontaler Versatz für 0.2 mm Normalwand
            lines.append((c0 + nh * off, c1 - c0))
        Q = np.array([line_intersect(lines[i - 1][0], lines[i - 1][1], lines[i][0], lines[i][1])
                      for i in range(4)])
        outer[j, :, :2], outer[j, :, 2] = C, z
        inner[j, :, :2], inner[j, :, 2] = Q, z
    return outer, inner


# ---------------- Mesh-Aufbau ----------------
tris = []  # Liste von (3,3)-Dreiecken, CCW von außen gesehen


def quad(a, b, c, d):
    tris.append((a, b, c))
    tris.append((a, c, d))


def ring_plate(z0, z1, twist, inner_holes, outer_holes, top_is_tube_side):
    """Ring mit 96 Öffnungen. tube-Seite: Loch = Röhren-Außenkontur; Gegenseite: Röhren-Innenkontur.
    Vertikale Lochwand = Innenkontur, verbunden mit Röhreninnenwand."""
    z_tube, z_free = (z1, z0) if top_is_tube_side else (z0, z1)
    up = 1 if top_is_tube_side else -1   # Normalenrichtung der Röhrenseite
    for k in range(N_SEG):
        th = k * PITCH + twist
        a, b = th - PITCH / 2, th + PITCH / 2
        S = np.array([[R_IN * np.cos(a), R_IN * np.sin(a)], [R_OUT * np.cos(a), R_OUT * np.sin(a)],
                      [R_OUT * np.cos(b), R_OUT * np.sin(b)], [R_IN * np.cos(b), R_IN * np.sin(b)]])
        O, I = outer_holes[k], inner_holes[k]
        P3 = lambda p, z: np.array([p[0], p[1], z])
        for i in range(4):
            j = (i + 1) % 4
            # Röhrenseite: Rahmen Sektor -> Außenkontur
            fa = [P3(S[i], z_tube), P3(S[j], z_tube), P3(O[j], z_tube), P3(O[i], z_tube)]
            # freie Seite: Rahmen Sektor -> Innenkontur
            fb = [P3(S[i], z_free), P3(S[j], z_free), P3(I[j], z_free), P3(I[i], z_free)]
            if up > 0:
                quad(*fa); quad(*fb[::-1])
            else:
                quad(*fa[::-1]); quad(*fb)
            # Lochwand (vertikal, Innenkontur) – Normale zeigt ins Loch
            w = [P3(I[i], z_free), P3(I[j], z_free), P3(I[j], z_tube), P3(I[i], z_tube)]
            quad(*(w[::-1] if up > 0 else w))
        # Außen- und Innenmantel des Rings
        for (r0, r1, sgn) in ((R_OUT, R_OUT, 1), (R_IN, R_IN, -1)):
            p0 = np.array([r0 * np.cos(a), r0 * np.sin(a)]); p1 = np.array([r1 * np.cos(b), r1 * np.sin(b)])
            m = [P3(p0, z0), P3(p1, z0), P3(p1, z1), P3(p0, z1)]
            quad(*(m if sgn > 0 else m[::-1]))


all_outer, all_inner = [], []
for k in range(N_SEG):
    o, i = tube_slices(k)
    all_outer.append(o); all_inner.append(i)
    for s in range(N_SLICES):
        for e in range(4):
            f = (e + 1) % 4
            quad(o[s, e], o[s, f], o[s + 1, f], o[s + 1, e])          # Außenwand
            quad(i[s, e], i[s + 1, e], i[s + 1, f], i[s, f])          # Innenwand

ring_plate(-PLATE_T, 0.0, 0.0,
           [x[0, :, :2] for x in all_inner], [x[0, :, :2] for x in all_outer], True)
ring_plate(HEIGHT, HEIGHT + PLATE_T, TWIST,
           [x[-1, :, :2] for x in all_inner], [x[-1, :, :2] for x in all_outer], False)

T = np.array(tris)  # (n,3,3)

# Orientierung aller Dreiecke prüfen/angleichen: Signiertes Volumen muss positiv sein
vol = np.einsum('ij,ij->i', T[:, 0], np.cross(T[:, 1], T[:, 2])).sum() / 6
if vol < 0:
    T = T[:, ::-1]
    vol = -vol

# Manifold-Check (Kanten verschweißen)
keys = np.round(T.reshape(-1, 3), 5)
uniq, idx = np.unique(keys, axis=0, return_inverse=True)
idx = idx.reshape(-1, 3)
edges = np.concatenate([idx[:, [0, 1]], idx[:, [1, 2]], idx[:, [2, 0]]])
directed = {tuple(e) for e in edges.tolist()}
bad = sum(1 for (a, b) in directed if (b, a) not in directed)
print(f"Dreiecke: {len(T):,}  Vertices: {len(uniq):,}  Volumen: {vol/1000:.2f} cm³  "
      f"offene/inkonsistente Kanten: {bad}  doppelte gerichtete Kanten: {len(edges)-len(directed)}")

# ---------------- STL (binär, mm) ----------------
def normals(T):
    n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    l = np.linalg.norm(n, axis=1, keepdims=True); l[l == 0] = 1
    return (n / l).astype(np.float32)


N = normals(T)
with open(OUT / "RingLampe.stl", "wb") as f:
    f.write(b"RingLampe AD110 ID100 H200 96 Roehren 0.2mm".ljust(80, b" "))
    f.write(struct.pack("<I", len(T)))
    rec = np.zeros(len(T), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = N, T.astype(np.float32)
    f.write(rec.tobytes())

# ---------------- GLB (Meter, Y-up, gelber transparenter Kunststoff) ----------------
pos = (T.reshape(-1, 3) * 0.001)[:, [0, 2, 1]] * np.array([1, 1, -1])  # Z-up mm -> Y-up m
pos = pos.astype(np.float32)
nrm = np.cross(pos[1::3] - pos[0::3], pos[2::3] - pos[0::3])
nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-20
nrm = np.repeat(nrm, 3, axis=0).astype(np.float32)
bin_pos, bin_nrm = pos.tobytes(), nrm.tobytes()
blob = bin_pos + bin_nrm
gltf = {
    "asset": {"version": "2.0", "generator": "generate_ringlampe.py"},
    "extensionsUsed": ["KHR_materials_transmission", "KHR_materials_ior", "KHR_materials_volume"],
    "scene": 0, "scenes": [{"nodes": [0]}],
    "nodes": [{"mesh": 0, "name": "RingLampe"}],
    "meshes": [{"name": "RingLampe", "primitives": [{"attributes": {"POSITION": 0, "NORMAL": 1}, "material": 0}]}],
    "materials": [{
        "name": "Kunststoff gelb transparent (PETG)",
        "doubleSided": True,
        "alphaMode": "BLEND",
        "pbrMetallicRoughness": {"baseColorFactor": [1.0, 0.82, 0.08, 0.55],
                                 "metallicFactor": 0.0, "roughnessFactor": 0.15},
        "emissiveFactor": [0.12, 0.09, 0.0],
        "extensions": {
            "KHR_materials_transmission": {"transmissionFactor": 0.85},
            "KHR_materials_ior": {"ior": 1.57},
            "KHR_materials_volume": {"thicknessFactor": 0.0002,
                                     "attenuationColor": [1.0, 0.8, 0.1], "attenuationDistance": 0.02},
        },
    }],
    "buffers": [{"byteLength": len(blob)}],
    "bufferViews": [
        {"buffer": 0, "byteOffset": 0, "byteLength": len(bin_pos), "target": 34962},
        {"buffer": 0, "byteOffset": len(bin_pos), "byteLength": len(bin_nrm), "target": 34962},
    ],
    "accessors": [
        {"bufferView": 0, "componentType": 5126, "count": len(pos), "type": "VEC3",
         "min": pos.min(0).tolist(), "max": pos.max(0).tolist()},
        {"bufferView": 1, "componentType": 5126, "count": len(nrm), "type": "VEC3"},
    ],
}
js = json.dumps(gltf, separators=(",", ":")).encode()
js += b" " * (-len(js) % 4)
blob += b"\0" * (-len(blob) % 4)
with open(OUT / "RingLampe.glb", "wb") as f:
    f.write(struct.pack("<III", 0x46546C67, 2, 12 + 8 + len(js) + 8 + len(blob)))
    f.write(struct.pack("<II", len(js), 0x4E4F534A)); f.write(js)
    f.write(struct.pack("<II", len(blob), 0x004E4942)); f.write(blob)

# ---------------- Technische Zeichnung (A3 quer) ----------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Polygon, Circle, Rectangle

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "lines.linewidth": 0.5})
LW_THICK, LW_THIN = 0.7, 0.25

waist_r_out = np.linalg.norm(all_outer[0][N_SLICES // 2, :, :2], axis=1).max()
waist_r_in = np.linalg.norm(all_outer[0][N_SLICES // 2, :, :2], axis=1).min()


def dim_h(ax, x0, x1, y, text, off=0, fs=7):
    ax.annotate("", (x0, y), (x1, y), arrowprops=dict(arrowstyle="<|-|>", lw=LW_THIN, color="k",
                                                     mutation_scale=6, shrinkA=0, shrinkB=0))
    ax.text((x0 + x1) / 2, y + off, text, ha="center", va="bottom", fontsize=fs)


def dim_v(ax, x, y0, y1, text, fs=7, side=-1):
    ax.annotate("", (x, y0), (x, y1), arrowprops=dict(arrowstyle="<|-|>", lw=LW_THIN, color="k",
                                                     mutation_scale=6, shrinkA=0, shrinkB=0))
    ax.text(x + side * 1.5, (y0 + y1) / 2, text, ha="right" if side < 0 else "left",
            va="center", rotation=90, fontsize=fs)


def setup(ax, title):
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=8, fontweight="bold", loc="left")


fig = plt.figure(figsize=(16.54, 11.69))  # A3 quer
fig.add_artist(Rectangle((0.015, 0.015), 0.97, 0.97, transform=fig.transFigure,
                         fill=False, lw=1.0))

# --- Vorderansicht ---
ax = fig.add_axes([0.03, 0.30, 0.33, 0.64])
setup(ax, "VORDERANSICHT (Röhrenkanten)")
for k in range(N_SEG):
    o = all_outer[k]
    for e in range(4):
        ax.plot(o[:, e, 0], o[:, e, 2], color="0.35", lw=0.15)
ax.add_patch(Rectangle((-R_OUT, -PLATE_T), 2 * R_OUT, PLATE_T, fill=False, lw=LW_THICK))
ax.add_patch(Rectangle((-R_OUT, HEIGHT), 2 * R_OUT, PLATE_T, fill=False, lw=LW_THICK))
ax.plot([0, 0], [-10, HEIGHT + 10], ls=(0, (12, 3, 2, 3)), color="k", lw=LW_THIN)
dim_v(ax, -R_OUT - 12, 0, HEIGHT, f"{HEIGHT:.0f}")
dim_v(ax, -R_OUT - 24, -PLATE_T, HEIGHT + PLATE_T, f"{HEIGHT + 2 * PLATE_T:.0f}")
dim_v(ax, R_OUT + 8, -PLATE_T, 0, f"{PLATE_T:.0f}", side=1)
dim_h(ax, -R_OUT, R_OUT, -14, f"Ø{2*R_OUT:.0f}")
dim_h(ax, -waist_r_out, waist_r_out, HEIGHT / 2 + 4, f"Ø{2*waist_r_out:.1f}", fs=6)
ax.text(R_OUT + 3, HEIGHT / 2, "A", fontsize=9, fontweight="bold", va="bottom")
ax.text(-R_OUT - 6, HEIGHT / 2, "A", fontsize=9, fontweight="bold", va="bottom")
ax.plot([-R_OUT - 4, -waist_r_out - 3], [HEIGHT / 2] * 2, color="k", lw=LW_THICK)
ax.plot([waist_r_out + 3, R_OUT + 4], [HEIGHT / 2] * 2, color="k", lw=LW_THICK)
ax.set_xlim(-R_OUT - 35, R_OUT + 20); ax.set_ylim(-25, HEIGHT + 12)


def ring_view(ax, twist, holes_o, holes_i, show_idx=True):
    ax.add_patch(Circle((0, 0), R_OUT, fill=False, lw=LW_THICK))
    ax.add_patch(Circle((0, 0), R_IN, fill=False, lw=LW_THICK))
    for k in range(N_SEG):
        ax.add_patch(Polygon(holes_o[k], closed=True, fill=False, lw=LW_THIN, color="0.4"))
        ax.add_patch(Polygon(holes_i[k], closed=True, fill=True, lw=LW_THIN, fc="0.85", ec="k"))
    for a in (twist, twist + np.pi / 2):
        ax.plot([-62 * np.cos(a), 62 * np.cos(a)], [-62 * np.sin(a), 62 * np.sin(a)],
                ls=(0, (12, 3, 2, 3)), color="k", lw=LW_THIN)
    if show_idx:
        p = seg_quad(0, twist).mean(0) * 1.18
        ax.text(*p, "LED 1", fontsize=6, ha="center", va="center")


# --- Draufsicht untere Basis ---
ax = fig.add_axes([0.36, 0.56, 0.30, 0.38])
setup(ax, "UNTERE BASIS – DRAUFSICHT")
ring_view(ax, 0.0, [x[0, :, :2] for x in all_outer], [x[0, :, :2] for x in all_inner])
dim_h(ax, -R_OUT, R_OUT, -64, f"Ø{2*R_OUT:.0f}")
dim_h(ax, -R_IN, R_IN, -5, f"Ø{2*R_IN:.0f}")
ax.text(0, 66, f"{N_SEG}× Segment, Teilung 360°/{N_SEG} = {360/N_SEG:.2f}°", ha="center", fontsize=6.5)
ax.add_patch(Circle(seg_quad(0).mean(0), 5, fill=False, lw=LW_THIN, ls="--"))
ax.text(R_OUT + 1, 4, "X", fontsize=9, fontweight="bold")
ax.set_xlim(-70, 70); ax.set_ylim(-70, 70)

# --- Draufsicht obere Basis ---
ax = fig.add_axes([0.66, 0.56, 0.30, 0.38])
setup(ax, "OBERE BASIS – DRAUFSICHT, um 120° verdreht")
ring_view(ax, TWIST, [x[-1, :, :2] for x in all_outer], [x[-1, :, :2] for x in all_inner])
ax.annotate("", seg_quad(0, TWIST).mean(0) * 0.8, seg_quad(0).mean(0) * 0.8,
            arrowprops=dict(arrowstyle="-|>", lw=LW_THIN, connectionstyle="arc3,rad=0.55", color="k"))
ax.text(28, 30, "120°", fontsize=7)
ax.plot(*(seg_quad(0).mean(0) * 0.8), "k+", ms=5)
ax.set_xlim(-70, 70); ax.set_ylim(-70, 70)

# --- Schnitt A-A ---
ax = fig.add_axes([0.36, 0.30, 0.20, 0.25])
setup(ax, "SCHNITT A–A  (z = 100)")
mid = N_SLICES // 2
for k in range(N_SEG):
    ax.add_patch(Polygon(all_outer[k][mid, :, :2], fill=True, fc="0.6", ec="k", lw=0.15))
    ax.add_patch(Polygon(all_inner[k][mid, :, :2], fill=True, fc="white", ec="k", lw=0.1))
dim_h(ax, -waist_r_out, waist_r_out, -33, f"Ø{2*waist_r_out:.2f}")
ax.text(0, 0, f"Taille: Ri {waist_r_in:.2f} / Ra {waist_r_out:.2f}\nSchnitt = Basis × 0.5, um 60° gedreht",
        ha="center", va="center", fontsize=5.5)
ax.set_xlim(-36, 36); ax.set_ylim(-38, 33)

# --- Detail X ---
ax = fig.add_axes([0.57, 0.27, 0.24, 0.28])
setup(ax, "DETAIL X (vergrößert) – Segment 1, Basis")
Rm = rot(-0.0)
O0, I0 = all_outer[0][0, :, :2], all_inner[0][0, :, :2]
O1 = all_outer[1][0, :, :2]
for poly, kw in ((O0, dict(fc="0.6", ec="k", lw=LW_THICK)), (I0, dict(fc="white", ec="k", lw=LW_THICK)),
                 (O1, dict(fc="0.6", ec="k", lw=LW_THIN)), (all_inner[1][0, :, :2], dict(fc="white", ec="k", lw=LW_THIN)),
                 (all_outer[-1][0, :, :2], dict(fc="0.6", ec="k", lw=LW_THIN)),
                 (all_inner[-1][0, :, :2], dict(fc="white", ec="k", lw=LW_THIN))):
    ax.add_patch(Polygon(poly, closed=True, **kw))
ax.plot([R_IN, R_IN], [-4, 4], "k", lw=LW_THIN); ax.plot([R_OUT, R_OUT], [-4, 4], "k", lw=LW_THIN)
dim_h(ax, R_IN, R_OUT, -4.6, "5.0 (Ring)", fs=6)
dim_h(ax, SEG_R_IN, SEG_R_OUT, -3.4, f"{SEG_R_OUT-SEG_R_IN:.1f}", fs=6)
wa = 2 * SEG_R_IN * np.sin(SEG_W / 2); wb = 2 * SEG_R_OUT * np.sin(SEG_W / 2)
dim_v(ax, SEG_R_IN - 0.25, O0[0, 1], O0[3, 1], f"{wa:.2f}", fs=6)
dim_v(ax, SEG_R_OUT + 0.25, O0[1, 1], O0[2, 1], f"{wb:.2f}", fs=6, side=1)
ax.annotate(f"Wand {WALL} (normal)", I0[1], (55.2, 4.0), fontsize=6,
            arrowprops=dict(arrowstyle="-", lw=LW_THIN))
ax.annotate(f"Spalt {GAP:.1f}", (SEG_R_IN + 0.15, (O0[3, 1] + O1[0, 1]) / 2), (47.0, 3.4), fontsize=6,
            arrowprops=dict(arrowstyle="-", lw=LW_THIN))
ax.text(52.5, 0, "LED\n2.0×2.0", ha="center", va="center", fontsize=5.5, color="0.3")
ax.set_xlim(46, 60); ax.set_ylim(-5.5, 5)

# --- Isometrie ---
ax3 = fig.add_axes([0.80, 0.27, 0.17, 0.29], projection="3d")
for k in range(0, N_SEG, 2):
    o = all_outer[k]
    ax3.plot(o[:, 1, 0], o[:, 1, 1], o[:, 1, 2], color=(0.85, 0.65, 0.0), lw=0.3)
for z, tw in ((0, 0), (HEIGHT, 0)):
    a = np.linspace(0, 2 * np.pi, 200)
    for r in (R_IN, R_OUT):
        ax3.plot(r * np.cos(a), r * np.sin(a), np.full_like(a, z), "k", lw=0.4)
ax3.set_box_aspect((1, 1, HEIGHT / 110)); ax3.set_axis_off(); ax3.view_init(18, -60)
ax3.set_title("ISOMETRIE (ohne Maßstab)", fontsize=8, fontweight="bold", loc="left")

# --- Hinweise & Schriftfeld ---
notes = [
    "HINWEISE",
    "1. Alle Maße in mm. Allgemeintoleranzen ISO 2768-m.",
    f"2. {N_SEG} Röhren, linear gelofted (Regelfläche) von Segment k (unten) zu Segment k (oben, +120°).",
    f"3. Röhrenwandstärke {WALL} mm normal zur Wandfläche; Röhren offen (Lichtkanal).",
    "4. Jeder Horizontalschnitt = Basis-Segmente × cos(φ/2), gedreht um φ/2 (φ = Verdrehwinkel z/H·120°)",
    "   → Röhren schneiden sich nicht; Mindestspalt in der Taille 0.3 mm.",
    f"5. LED-Öffnung (lichte Weite Basis) ca. {np.ptp(all_inner[0][0,:,:2]@rot(-0)[:,:].T,axis=0).min():.1f} × "
    f"{np.linalg.norm(all_inner[0][0,1]-all_inner[0][0,0]):.1f} mm, geeignet für LEDs 2020 (2.0×2.0 mm).",
    "6. Material: PETG / PMMA transluzent gelb. Fertigung: SLA/Resin empfohlen (0.2-mm-Wände).",
    f"7. Volumen ca. {vol/1000:.1f} cm³.",
]
fig.text(0.03, 0.235, "\n".join(notes), fontsize=7, va="top", family="DejaVu Sans", linespacing=1.6)

tb = fig.add_axes([0.60, 0.025, 0.375, 0.16]); tb.axis("off")
tb.set_xlim(0, 100); tb.set_ylim(0, 40)
for y in (0, 10, 20, 30, 40):
    tb.plot([0, 100], [y, y], "k", lw=0.6)
for x in (0, 100):
    tb.plot([x, x], [0, 40], "k", lw=0.6)
tb.plot([50, 50], [0, 30], "k", lw=0.6)
rows = [
    (2, 35, "Benennung:", "RingLampe – 96-LED Röhrenkörper (120° verdreht)", 9),
    (2, 25, "Werkstoff:", "Kunststoff gelb transparent (PETG/PMMA)", 7),
    (52, 25, "Maßstab:", "unmaßstäblich – Maße gelten, A3", 7),
    (2, 15, "Datum:", date.today().isoformat(), 7),
    (52, 15, "Zeichnungs-Nr.:", "RL-001  Blatt 1/1", 7),
    (2, 5, "Dateien:", "RingLampe.stl / RingLampe.glb", 7),
    (52, 5, "Projektion:", "Methode 1 (ISO E)", 7),
]
for x, y, lab, val, fs in rows:
    tb.text(x, y + 2.2, lab, fontsize=5.5, color="0.35")
    tb.text(x, y - 2.5, val, fontsize=fs, fontweight="bold" if fs > 7 else "normal")

with PdfPages(OUT / "RingLampe_Zeichnung.pdf") as pdf:
    pdf.savefig(fig)

# Kontrolle: minimaler Spalt zwischen Nachbarröhren (Außenkonturen) über alle Höhen
gaps = []
for j in range(N_SLICES + 1):
    A, B = all_outer[0][j, :, :2], all_outer[1][j, :, :2]
    # Abstand Kante A(2->3) zu Kante B(0->1) über Eckpunkte (Kanten nahezu parallel)
    def pdist(p, q0, q1):
        d = q1 - q0; t = np.clip(np.dot(p - q0, d) / np.dot(d, d), 0, 1); return np.linalg.norm(p - q0 - t * d)
    gaps.append(min(pdist(A[2], B[0], B[1]), pdist(A[3], B[0], B[1]), pdist(B[0], A[3], A[2]), pdist(B[1], A[3], A[2])))
print(f"Min. Spalt zwischen Nachbarröhren: {min(gaps):.3f} mm (Basis {gaps[0]:.3f} mm)")
print(f"Taille: Ø außen {2*waist_r_out:.2f} mm")
