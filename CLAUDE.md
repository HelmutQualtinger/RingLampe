# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

RingLampe: a procedurally generated 3D lamp body. Two ring bases (OD 110 / ID 100 mm, 2 mm thick), each with 96 rectangular LED openings, sit 200 mm apart. The top ring is rotated 120°. 96 tubes with 0.2 mm walls are lofted from each bottom segment to its matching top segment. The original German spec is in `prompt.md`. Code comments, UI text and drawing labels are in German.

## Commands

```bash
python3 generate_ringlampe.py   # regenerates RingLampe.stl, RingLampe.glb, RingLampe_Zeichnung.pdf
```

- Dependencies: `numpy`, `matplotlib`. There is no build system, test suite or linter.
- The script runs top to bottom (no `main`, no CLI args). It always writes the outputs next to itself.
- Use the printed output as the check: the triangle count and volume, `offene/inkonsistente Kanten` (must be 0 for a watertight mesh) and `Min. Spalt zwischen Nachbarröhren` (must be > 0, so tubes don't intersect).
- Viewer: `ringlampe_viewer.html` loads three.js 0.160 as ES modules from jsDelivr, so serve it over HTTP: `python3 -m http.server 8000`, then open `http://localhost:8000/ringlampe_viewer.html`. Don't use port 5000.

## Architecture

**`generate_ringlampe.py`** builds one triangle soup `tris` in mm, Z-up:
- `seg_quad(k, twist)` gives the outer footprint of segment k. `tube_slices(k)` linearly interpolates bottom quad P to the rotated quad RP across `N_SLICES` z-levels. It returns the outer contour and an inner contour. The inner contour is offset by `WALL / h` per edge, so the wall is 0.2 mm *normal to the ruled surface*, not horizontally.
- Why tubes can't intersect: every horizontal cross-section is a similarity transform (rotation + scaling by cos(φ/2)) of the base layout. So disjoint base segments stay disjoint. Changing the interpolation from linear breaks this guarantee.
- `ring_plate(...)` stitches each plate to the tube ends. On the tube side, the holes match the tube's outer contour. On the free side, they match the inner contour. The vertical hole walls connect to the tube's inner wall. This shared-vertex stitching is what keeps the mesh manifold, so keep tube end contours and plate hole contours identical.
- After assembly, the script flips the orientation if the signed volume is negative, then welds vertices (rounded to 1e-5) to count open edges.
- Exports:
  - **STL**: binary, in mm.
  - **GLB**: hand-written glTF 2.0, converted to meters and Y-up (`[x, z, -y]`). Its material is yellow transparent PETG via `KHR_materials_transmission`, `KHR_materials_ior` and `KHR_materials_volume`.
  - **PDF**: an A3 technical drawing made with matplotlib. It has front, top-bottom and top-top views, section A–A at z=100, detail X, an isometric view, notes and a title block. Some note text is hard-coded (for example "Mindestspalt 0.3 mm" and the LED size), so update it when you change the parameters.

**`ringlampe_viewer.html`** is a standalone three.js viewer. It rebuilds the geometry in JS and does **not** load the GLB/STL. Its parameter block (`R_IN`, `R_OUT`, `N`, `H`, `TWIST`, `SEG_R_*`, `GAP`, …) duplicates the Python constants, so keep the two in sync by hand.
- The tubes use a custom `ShaderMaterial`. Per-vertex `aSeg` and `aT` attributes drive a rotating rainbow that travels upward with a delay and fades with height. RGB LEDs are an `InstancedMesh` under the bottom plate, and 8 point lights give colored ambient fill. Bloom comes from `EffectComposer` + `UnrealBloomPass`. A "tumble" mode rotates the model pivot about 3 axes.
- `?record` mode renders deterministically at 1080×1350 (4:5), 30 fps, for 60 s, with a text overlay and fades. `window.startRecording(endpoint='/frame/')` POSTs each JPEG frame to a local server, and `window.__rec` reports progress. That frame-receiving server isn't in this repo, so you have to provide it when recording.

Generated artifacts (`RingLampe.stl`, `RingLampe.glb`, `RingLampe_Zeichnung.pdf`) are committed outputs. Regenerate them rather than editing them.
