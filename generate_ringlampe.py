#!/usr/bin/env python3
"""RingLampe – prozedurale Erzeugung (STL, GLB, technische Zeichnung PDF).

Geometrie (Einheiten mm):
  * 2 ringförmige Basen, AD 110 / ID 100, Dicke PLATE_T, mit je 48 runden LED-Öffnungen
  * obere Basis 200 mm über der unteren, um 120° verdreht
  * 48 runde Röhren, linear gelofted von der unteren zur oberen Öffnung, Wandstärke 0.2 mm (normal zur Wand)
  * Jeder Horizontalschnitt des Röhrenbündels ist eine Ähnlichkeitsabbildung (Drehung+Skalierung)
    der Basis-Kreise -> Querschnitte bleiben Kreise, Röhren schneiden sich nie, solange die
    Basis-Kreise disjunkt sind.
"""
import json
import struct
from datetime import date
from pathlib import Path

import numpy as np

OUT = Path(__file__).parent

# ---------------- Parameter ----------------
R_IN, R_OUT = 50.0, 55.0          # Ring ID 100 / AD 110
N_SEG = 48                        # LEDs / Segmente
HEIGHT = 200.0                    # Abstand der Basen
TWIST = np.radians(120.0)         # Verdrehung obere Basis
WALL = 0.2                        # Röhrenwandstärke
PLATE_T = 2.0                     # Dicke der Basisringe
TUBE_RC = 52.5                    # Teilkreisradius der Röhrenmitten
TUBE_D = 2.8                      # Röhren-Außendurchmesser an der Basis
N_CIRC = 16                       # Polygonpunkte je Röhrenquerschnitt
N_SLICES = 40                     # Unterteilung der Röhren in z

PITCH = 2 * np.pi / N_SEG
GAP = 2 * TUBE_RC * np.sin(PITCH / 2) - TUBE_D   # Spalt zwischen Nachbarröhren an der Basis


def rot(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s], [s, c]])


def tube_center(k, twist=0.0):
    th = k * PITCH + twist
    return TUBE_RC * np.array([np.cos(th), np.sin(th)])


def tube_circle(k, twist=0.0):
    """Außenkontur (Footprint) von Röhre k, CCW, N_CIRC Punkte. Startwinkel dreht mit,
    damit tube_circle(k, φ) = rot(φ) · tube_circle(k, 0) punktweise gilt."""
    th = k * PITCH + twist
    psi = th + 2 * np.pi * np.arange(N_CIRC) / N_CIRC
    return tube_center(k, twist) + TUBE_D / 2 * np.c_[np.cos(psi), np.sin(psi)]


def tube_slices(k):
    """Liefert outer[N+1,M,3], inner[N+1,M,3] für Röhre k."""
    P = tube_circle(k)
    RP = P @ rot(TWIST).T
    D = np.c_[RP - P, np.full(N_CIRC, HEIGHT)]          # Erzeugende je Konturpunkt
    ts = np.linspace(0, 1, N_SLICES + 1)
    outer = np.zeros((len(ts), N_CIRC, 3))
    inner = np.zeros((len(ts), N_CIRC, 3))
    for j, t in enumerate(ts):
        C = (1 - t) * P + t * RP
        cen = C.mean(axis=0)
        e = np.c_[np.roll(C, -1, 0) - np.roll(C, 1, 0), np.zeros(N_CIRC)]   # Tangente
        n = np.cross(e, D)
        n /= np.linalg.norm(n, axis=1, keepdims=True)
        h = np.hypot(n[:, 0], n[:, 1])
        nh = n[:, :2] / h[:, None]
        nh *= np.sign(np.einsum('ij,ij->i', cen - C, nh))[:, None]           # nach innen
        Q = C + nh * (WALL / h)[:, None]       # horizontaler Versatz für 0.2 mm Normalwand
        outer[j, :, :2], outer[j, :, 2] = C, HEIGHT * t
        inner[j, :, :2], inner[j, :, 2] = Q, HEIGHT * t
    return outer, inner


# ---------------- Mesh-Aufbau ----------------
tris = []  # Liste von (3,3)-Dreiecken, CCW von außen gesehen


def quad(a, b, c, d):
    tris.append((a, b, c))
    tris.append((a, c, d))


def zipper(A, B):
    """Trianguliert den Ring zwischen äußerer Schleife A und innerer Schleife B (beide CCW,
    sternförmig um B-Mitte). Dreiecke CCW von +z gesehen, nur Originalpunkte (keine T-Stöße)."""
    c = B.mean(0)
    ref = np.arctan2(*(A[0] - c)[::-1])
    ang = lambda P: np.mod(np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0]) - ref, 2 * np.pi)
    B = np.roll(B, -np.argmin(ang(B)), 0)
    aa, ab = np.append(ang(A), 2 * np.pi), ang(B)
    ab = np.append(ab, ab[0] + 2 * np.pi)
    A2, B2 = np.vstack([A, A[:1]]), np.vstack([B, B[:1]])
    out, i, j = [], 0, 0
    while i < len(A) or j < len(B):
        if j == len(B) or (i < len(A) and aa[i + 1] <= ab[j + 1]):
            out.append((A2[i], A2[i + 1], B2[j])); i += 1
        else:
            out.append((A2[i], B2[j + 1], B2[j])); j += 1
    return out


def ring_plate(z0, z1, twist, inner_holes, outer_holes, top_is_tube_side):
    """Ring mit N_SEG Öffnungen. tube-Seite: Loch = Röhren-Außenkontur; Gegenseite: Röhren-Innenkontur.
    Vertikale Lochwand = Innenkontur, verbunden mit Röhreninnenwand."""
    z_tube, z_free = (z1, z0) if top_is_tube_side else (z0, z1)
    up = 1 if top_is_tube_side else -1   # Normalenrichtung der Röhrenseite
    P3 = lambda p, z: np.array([p[0], p[1], z])
    for k in range(N_SEG):
        th = k * PITCH + twist
        a, b = th - PITCH / 2, th + PITCH / 2
        S = np.array([[R_IN * np.cos(a), R_IN * np.sin(a)], [R_OUT * np.cos(a), R_OUT * np.sin(a)],
                      [R_OUT * np.cos(b), R_OUT * np.sin(b)], [R_IN * np.cos(b), R_IN * np.sin(b)]])
        O, I = outer_holes[k], inner_holes[k]
        # Röhrenseite: Rahmen Sektor -> Außenkontur; freie Seite: Rahmen Sektor -> Innenkontur
        for loop, z, sgn in ((O, z_tube, up), (I, z_free, -up)):
            for t in zipper(S, loop):
                f = [P3(p, z) for p in t]
                tris.append(tuple(f if sgn > 0 else f[::-1]))
        for i in range(N_CIRC):
            j = (i + 1) % N_CIRC
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
        for e in range(N_CIRC):
            f = (e + 1) % N_CIRC
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


# Kontrolle: minimaler Spalt zwischen Nachbarröhren (Außenkonturen) über alle Höhen
def poly_dist(A, B):
    """Minimaler Abstand zweier geschlossener Polygone (Punkt-zu-Kante, beide Richtungen)."""
    def pe(P, Q):
        q0, d = Q, np.roll(Q, -1, 0) - Q
        w = P[:, None] - q0[None]
        t = np.clip(np.einsum('pqi,qi->pq', w, d) / np.einsum('qi,qi->q', d, d), 0, 1)
        return np.linalg.norm(w - t[..., None] * d[None], axis=2).min()
    return min(pe(A, B), pe(B, A))


gaps = [poly_dist(all_outer[0][j, :, :2], all_outer[1][j, :, :2]) for j in range(N_SLICES + 1)]
min_gap = min(gaps)

# ---------------- STL (binär, mm) ----------------
def normals(T):
    n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    l = np.linalg.norm(n, axis=1, keepdims=True); l[l == 0] = 1
    return (n / l).astype(np.float32)


N = normals(T)
with open(OUT / "RingLampe.stl", "wb") as f:
    f.write(f"RingLampe AD110 ID100 H200 {N_SEG} runde Roehren 0.2mm".encode().ljust(80, b" "))
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
    for e in range(0, N_CIRC, N_CIRC // 4):
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
        p = tube_center(0, twist) * 1.18
        ax.text(*p, "LED 1", fontsize=6, ha="center", va="center")


# --- Draufsicht untere Basis ---
ax = fig.add_axes([0.36, 0.56, 0.30, 0.38])
setup(ax, "UNTERE BASIS – DRAUFSICHT")
ring_view(ax, 0.0, [x[0, :, :2] for x in all_outer], [x[0, :, :2] for x in all_inner])
dim_h(ax, -R_OUT, R_OUT, -64, f"Ø{2*R_OUT:.0f}")
dim_h(ax, -R_IN, R_IN, -5, f"Ø{2*R_IN:.0f}")
ax.text(0, 66, f"{N_SEG}× Bohrung Ø{TUBE_D:.1f} auf Teilkreis Ø{2*TUBE_RC:.0f}, "
        f"Teilung 360°/{N_SEG} = {360/N_SEG:.2f}°", ha="center", fontsize=6.5)
ax.add_patch(Circle(tube_center(0), 5, fill=False, lw=LW_THIN, ls="--"))
ax.text(R_OUT + 1, 4, "X", fontsize=9, fontweight="bold")
ax.set_xlim(-70, 70); ax.set_ylim(-70, 70)

# --- Draufsicht obere Basis ---
ax = fig.add_axes([0.66, 0.56, 0.30, 0.38])
setup(ax, "OBERE BASIS – DRAUFSICHT, um 120° verdreht")
ring_view(ax, TWIST, [x[-1, :, :2] for x in all_outer], [x[-1, :, :2] for x in all_inner])
ax.annotate("", tube_center(0, TWIST) * 0.8, tube_center(0) * 0.8,
            arrowprops=dict(arrowstyle="-|>", lw=LW_THIN, connectionstyle="arc3,rad=0.55", color="k"))
ax.text(28, 30, "120°", fontsize=7)
ax.plot(*(tube_center(0) * 0.8), "k+", ms=5)
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
YL = TUBE_RC * np.sin(PITCH) + TUBE_D / 2 + 0.4          # Nachbarröhren noch sichtbar
ax.plot([R_IN, R_IN], [-YL, YL], "k", lw=LW_THIN); ax.plot([R_OUT, R_OUT], [-YL, YL], "k", lw=LW_THIN)
dim_h(ax, R_IN, R_OUT, -TUBE_D / 2 - 1.6, "5.0 (Ring)", fs=6)
c0 = tube_center(0)
ax.annotate(f"Ø{TUBE_D:.1f} außen", c0 + TUBE_D / 2 * np.array([0.71, -0.71]), (55.6, -1.8), fontsize=6,
            arrowprops=dict(arrowstyle="-", lw=LW_THIN))
d_in = 2 * np.linalg.norm(I0 - c0, axis=1).min()
ax.plot([c0[0] - 2.2, c0[0] + 2.2], [0, 0], ls=(0, (8, 2, 1, 2)), color="k", lw=LW_THIN)
ax.add_patch(Rectangle((c0[0] - 1, -1), 2, 2, fill=False, ls="--", lw=LW_THIN, ec="0.45"))
ax.annotate(f"Wand {WALL} (normal)", O0[N_CIRC // 8], (55.2, 4.0), fontsize=6,
            arrowprops=dict(arrowstyle="-", lw=LW_THIN))
ax.annotate(f"Spalt {GAP:.2f}", (TUBE_RC, (O0[:, 1].max() + O1[:, 1].min()) / 2), (46.3, (O0[:, 1].max() + O1[:, 1].min()) / 2 + 0.8), fontsize=6,
            arrowprops=dict(arrowstyle="-", lw=LW_THIN))
ax.annotate(f"lichte Weite Ø{d_in:.1f}", (c0[0] - d_in / 2 * 0.7, -d_in / 2 * 0.7), (46.3, -2.4), fontsize=6,
            arrowprops=dict(arrowstyle="-", lw=LW_THIN))
ax.annotate("LED 2.0×2.0 (unter Basis)", (c0[0] + 1, 1), (55.8, 1.6), fontsize=5.5, color="0.3",
            arrowprops=dict(arrowstyle="-", lw=LW_THIN, color="0.45"))
ax.text(46.3, 0.25, f"Teilkreis R{TUBE_RC:.1f}", ha="left", va="bottom", fontsize=5.5)
ax.plot([48.6, c0[0] - 2.2], [0, 0], color="k", lw=LW_THIN)
ax.set_xlim(46, 60); ax.set_ylim(-YL, YL)

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
    f"2. {N_SEG} runde Röhren Ø{TUBE_D:.1f} (Basis), linear gelofted (Regelfläche) von Bohrung k (unten)"
    " zu Bohrung k (oben, +120°).",
    f"3. Röhrenwandstärke {WALL} mm normal zur Wandfläche; Röhren offen (Lichtkanal).",
    "4. Jeder Horizontalschnitt = Basis-Kreise × cos(φ/2), gedreht um φ/2 (φ = Verdrehwinkel z/H·120°)",
    f"   → Querschnitte bleiben Kreise, Röhren schneiden sich nicht; Mindestspalt in der Taille {min_gap:.2f} mm.",
    f"5. LED-Öffnung rund, lichte Weite Basis ca. Ø{d_in:.1f} mm, über LEDs 2020 (2.0×2.0 mm).",
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
    (2, 35, "Benennung:", f"RingLampe – {N_SEG}-LED Rundröhrenkörper (120° verdreht)", 9),
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

print(f"Min. Spalt zwischen Nachbarröhren: {min_gap:.3f} mm (Basis {gaps[0]:.3f} mm)")
print(f"Lichte Weite Basis: Ø{d_in:.2f} mm")
print(f"Taille: Ø außen {2*waist_r_out:.2f} mm")
