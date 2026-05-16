#!/usr/bin/env python3
"""Create a PyMOL-like browser viewer from an H5 cache plus source structures.

The H5 file provides voxel labels and APBS grids. The protein and ligand source
files provide atom/bond level context for a molecule-style view.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import h5py
import numpy as np


DEFAULT_DATA_ROOT = Path("/Users/tevfik/Sandbox/github/PHD/data")
APBS_CANDIDATES = [
    "electrostatic_grid_v1_ligand_proximal_chains_7A_raw",
    "electrostatic_grid",
    "electrostatic_grid_v2_full_protein_raw",
]
PROVIDED_LABEL_CANDIDATES = [
    "binding_site_in_dataset",
    "binding_site_fpocket_selected",
    "binding_site_cavity6",
    "binding_site_cavityALL",
    "binding_site_site",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a PyMOL-like molecular HTML viewer.")
    parser.add_argument("h5_path", help="H5 cache file.")
    parser.add_argument("--protein", default=None, help="Protein structure file, preferably PDB.")
    parser.add_argument("--ligand", default=None, help="Ligand structure file, PDB/MOL2/SDF.")
    parser.add_argument("--data-root", default=str(DEFAULT_DATA_ROOT), help="Root used for automatic source lookup.")
    parser.add_argument("--output", default=None, help="Output HTML path.")
    parser.add_argument("--apbs-feature", default=None, help="APBS feature to visualize. Defaults to first available raw APBS.")
    parser.add_argument("--label", default=None, help="Dataset/provided label name. Defaults to first available provided label.")
    parser.add_argument("--max-apbs-points", type=int, default=1800, help="Maximum APBS positive and negative points each.")
    parser.add_argument("--max-label-points", type=int, default=2000, help="Maximum label points.")
    parser.add_argument("--apbs-percentile", type=float, default=98.5, help="Absolute APBS percentile used for point clouds.")
    parser.add_argument("--sphere-radius", type=float, default=0.7, help="Sphere radius for voxel-derived points.")
    return parser.parse_args()


def case_name_from_h5(path: Path) -> str:
    stem = path.stem
    for token in ("_box36", "_box72", "_box161", "_span70", "_span120", "_span160"):
        if token in stem:
            stem = stem.split(token)[0]
    return stem


def first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def auto_source_paths(case_name: str, data_root: Path) -> tuple[Path | None, Path | None]:
    candidates = [
        (
            data_root / "external_benchmarks/puresnet_prepared/bu48_puresnet" / case_name / "protein.pdb",
            data_root / "external_benchmarks/puresnet_prepared/bu48_puresnet" / case_name / "ligand.pdb",
        ),
        (
            data_root / "external_benchmarks/puresnet_prepared/coach420_puresnet" / case_name / "protein.pdb",
            data_root / "external_benchmarks/puresnet_prepared/coach420_puresnet" / case_name / "ligand.pdb",
        ),
        (
            data_root / "pdbbind/refined-set" / case_name / f"{case_name}_protein.pdb",
            data_root / "pdbbind/refined-set" / case_name / f"{case_name}_ligand.mol2",
        ),
        (
            data_root / "data_compressed/refined-set" / case_name / f"{case_name}_protein.pdb",
            data_root / "data_compressed/refined-set" / case_name / f"{case_name}_ligand.mol2",
        ),
        (
            data_root / "scPDB_converted" / case_name / f"{case_name}_protein.pdb",
            data_root / "scPDB_converted" / case_name / f"{case_name}_ligand.mol2",
        ),
        (
            data_root / "scPDB" / case_name / "protein.mol2",
            data_root / "scPDB" / case_name / "ligand.mol2",
        ),
    ]
    for protein, ligand in candidates:
        if protein.exists() and ligand.exists():
            return protein, ligand
    protein = first_existing([protein for protein, _ in candidates])
    ligand = first_existing([ligand for _, ligand in candidates])
    return protein, ligand


def read_text(path: Path | None) -> str | None:
    if path is None:
        return None
    return path.read_text(errors="replace")


def structure_format(path: Path | None) -> str:
    if path is None:
        return "pdb"
    suffix = path.suffix.lower().lstrip(".")
    if suffix == "ent":
        return "pdb"
    if suffix in {"mol2", "sdf", "pdb"}:
        return suffix
    return "pdb"


def choose_dataset(h5f: h5py.File, group: str, names: list[str], explicit: str | None = None) -> tuple[str | None, np.ndarray | None]:
    if explicit:
        path = f"{group}/{explicit}"
        if path in h5f:
            return explicit, h5f[path][:]
        return explicit, None
    for name in names:
        path = f"{group}/{name}"
        if path in h5f:
            return name, h5f[path][:]
    return None, None


def grid_to_world(coords: np.ndarray, h5f: h5py.File) -> np.ndarray:
    origin = np.asarray(h5f.attrs.get("grid_origin", h5f.attrs.get("electrostatic_grid_min", [0, 0, 0])), dtype=np.float32)
    resolution = float(h5f.attrs.get("resolution", 1.0))
    return origin + coords.astype(np.float32) * resolution


def sample_coords(coords: np.ndarray, values: np.ndarray | None, max_points: int) -> tuple[np.ndarray, np.ndarray | None]:
    if len(coords) <= max_points:
        return coords, values
    step = math.ceil(len(coords) / max_points)
    coords = coords[::step]
    if values is not None:
        values = values[::step]
    return coords, values


def apbs_point_clouds(apbs: np.ndarray, h5f: h5py.File, percentile: float, max_points: int):
    finite = apbs[np.isfinite(apbs)]
    if finite.size == 0:
        return [], []
    threshold = float(np.percentile(np.abs(finite), percentile))
    if threshold <= 0:
        threshold = float(np.max(np.abs(finite)))
    if threshold <= 0:
        return [], []

    pos_coords = np.argwhere(apbs >= threshold)
    neg_coords = np.argwhere(apbs <= -threshold)
    pos_vals = apbs[tuple(pos_coords.T)] if len(pos_coords) else np.array([])
    neg_vals = apbs[tuple(neg_coords.T)] if len(neg_coords) else np.array([])
    pos_coords, pos_vals = sample_coords(pos_coords, pos_vals, max_points)
    neg_coords, neg_vals = sample_coords(neg_coords, neg_vals, max_points)
    pos_world = grid_to_world(pos_coords, h5f) if len(pos_coords) else np.empty((0, 3), dtype=np.float32)
    neg_world = grid_to_world(neg_coords, h5f) if len(neg_coords) else np.empty((0, 3), dtype=np.float32)
    return points_payload(pos_world, pos_vals), points_payload(neg_world, neg_vals)


def mask_point_cloud(mask: np.ndarray | None, h5f: h5py.File, max_points: int):
    if mask is None:
        return []
    coords = np.argwhere(mask > 0)
    coords, _ = sample_coords(coords, None, max_points)
    world = grid_to_world(coords, h5f) if len(coords) else np.empty((0, 3), dtype=np.float32)
    return points_payload(world, None)


def points_payload(coords: np.ndarray, values: np.ndarray | None):
    out = []
    for idx, point in enumerate(coords):
        item = [round(float(point[0]), 3), round(float(point[1]), 3), round(float(point[2]), 3)]
        if values is not None and len(values):
            item.append(round(float(values[idx]), 3))
        out.append(item)
    return out


def js_json(value) -> str:
    return json.dumps(value, ensure_ascii=False)


def build_html(payload: dict) -> str:
    return f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{payload["title"]}</title>
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <style>
    body {{ margin:0; font-family:Arial, sans-serif; background:#0b1020; color:#e5e7eb; }}
    header {{ padding:14px 18px; background:#111827; border-bottom:1px solid #334155; }}
    h1 {{ margin:0; font-size:20px; }}
    .meta {{ color:#94a3b8; font-size:12px; margin-top:6px; }}
    #wrap {{ display:grid; grid-template-columns:280px 1fr; height:calc(100vh - 67px); }}
    #controls {{ padding:16px; background:#0f172a; border-right:1px solid #334155; overflow:auto; }}
    #viewer {{ position:relative; width:100%; height:100%; }}
    label {{ display:flex; gap:8px; align-items:center; margin:9px 0; font-size:13px; }}
    button {{ width:100%; margin:7px 0; padding:8px 10px; border:1px solid #475569; background:#1e293b; color:#e5e7eb; border-radius:6px; cursor:pointer; }}
    .legend {{ margin-top:14px; font-size:12px; color:#cbd5e1; line-height:1.55; }}
    .swatch {{ display:inline-block; width:10px; height:10px; border-radius:50%; margin-right:6px; }}
    .warn {{ margin-top:12px; padding:10px; background:#451a03; border:1px solid #92400e; border-radius:6px; color:#fed7aa; font-size:12px; }}
  </style>
</head>
<body>
  <header>
    <h1>{payload["title"]}</h1>
    <div class="meta">{payload["meta"]}</div>
  </header>
  <div id="wrap">
    <aside id="controls">
      <button onclick="resetView()">Reset view</button>
      <button onclick="setProteinStyle('cartoon')">Protein cartoon</button>
      <button onclick="setProteinStyle('surface')">Protein surface</button>
      <button onclick="setProteinStyle('stick')">Protein sticks</button>
      <label><input type="checkbox" id="ligandToggle" checked onchange="redraw()"> Ligand</label>
      <label><input type="checkbox" id="surfaceToggle" checked onchange="redraw()"> Transparent protein surface</label>
      <label><input type="checkbox" id="apbsPosToggle" checked onchange="redraw()"> Positive APBS cloud</label>
      <label><input type="checkbox" id="apbsNegToggle" checked onchange="redraw()"> Negative APBS cloud</label>
      <label><input type="checkbox" id="calcToggle" checked onchange="redraw()"> Calculated pocket label</label>
      <label><input type="checkbox" id="providedToggle" checked onchange="redraw()"> Dataset pocket label</label>
      <div class="legend">
        <div><span class="swatch" style="background:#2563eb"></span>Positive APBS</div>
        <div><span class="swatch" style="background:#dc2626"></span>Negative APBS</div>
        <div><span class="swatch" style="background:#22c55e"></span>Ligand</div>
        <div><span class="swatch" style="background:#f97316"></span>Calculated label</div>
        <div><span class="swatch" style="background:#a855f7"></span>Dataset label</div>
      </div>
      {payload["warning_html"]}
    </aside>
    <main id="viewer"></main>
  </div>
  <script>
    const payload = {js_json(payload)};
    const viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "#0b1020" }});
    let proteinModel = null;
    let ligandModel = null;
    let proteinStyle = "cartoon";

    function addSphereCloud(points, color, radius, alpha) {{
      for (const p of points) {{
        viewer.addSphere({{ center: {{x:p[0], y:p[1], z:p[2]}}, radius, color, alpha }});
      }}
    }}

    function applyProteinStyle() {{
      if (!proteinModel) return;
      if (proteinStyle === "cartoon") {{
        viewer.setStyle({{model: proteinModel}}, {{cartoon: {{color: "spectrum"}}}});
      }} else if (proteinStyle === "surface") {{
        viewer.setStyle({{model: proteinModel}}, {{}});
        viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity:0.72, color:"#dbeafe"}}, {{model: proteinModel}});
      }} else {{
        viewer.setStyle({{model: proteinModel}}, {{stick: {{radius:0.14, colorscheme:"Jmol"}}}});
      }}
    }}

    function setProteinStyle(style) {{
      proteinStyle = style;
      redraw();
    }}

    function redraw() {{
      viewer.clear();
      if (payload.proteinData) {{
        proteinModel = viewer.addModel(payload.proteinData, payload.proteinFormat);
        applyProteinStyle();
        if (document.getElementById("surfaceToggle").checked && proteinStyle !== "surface") {{
          viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity:0.18, color:"#e2e8f0"}}, {{model: proteinModel}});
        }}
      }}
      if (payload.ligandData && document.getElementById("ligandToggle").checked) {{
        ligandModel = viewer.addModel(payload.ligandData, payload.ligandFormat);
        viewer.setStyle({{model: ligandModel}}, {{stick: {{radius:0.26, colorscheme:"greenCarbon"}}}});
      }}
      if (document.getElementById("apbsPosToggle").checked) {{
        addSphereCloud(payload.apbsPositive, "#2563eb", payload.sphereRadius, 0.42);
      }}
      if (document.getElementById("apbsNegToggle").checked) {{
        addSphereCloud(payload.apbsNegative, "#dc2626", payload.sphereRadius, 0.42);
      }}
      if (document.getElementById("calcToggle").checked) {{
        addSphereCloud(payload.calculatedLabel, "#f97316", payload.sphereRadius * 0.85, 0.55);
      }}
      if (document.getElementById("providedToggle").checked) {{
        addSphereCloud(payload.providedLabel, "#a855f7", payload.sphereRadius * 0.85, 0.55);
      }}
      viewer.zoomTo();
      viewer.render();
    }}

    function resetView() {{
      viewer.zoomTo();
      viewer.render();
    }}

    redraw();
  </script>
</body>
</html>
"""


def main() -> None:
    args = parse_args()
    h5_path = Path(args.h5_path).expanduser().resolve()
    data_root = Path(args.data_root).expanduser().resolve()
    case_name = case_name_from_h5(h5_path)

    protein_path = Path(args.protein).expanduser().resolve() if args.protein else None
    ligand_path = Path(args.ligand).expanduser().resolve() if args.ligand else None
    if protein_path is None or ligand_path is None:
        auto_protein, auto_ligand = auto_source_paths(case_name, data_root)
        protein_path = protein_path or auto_protein
        ligand_path = ligand_path or auto_ligand

    with h5py.File(h5_path, "r") as h5f:
        apbs_name, apbs = choose_dataset(h5f, "features", APBS_CANDIDATES, args.apbs_feature)
        if apbs is None and "auxiliary" in h5f:
            apbs_name, apbs = choose_dataset(h5f, "auxiliary", APBS_CANDIDATES, args.apbs_feature)
        provided_name, provided = choose_dataset(h5f, "label", PROVIDED_LABEL_CANDIDATES, args.label)
        calculated = h5f["label/binding_site_calculated"][:] if "label/binding_site_calculated" in h5f else None
        apbs_positive, apbs_negative = apbs_point_clouds(
            apbs,
            h5f,
            percentile=args.apbs_percentile,
            max_points=args.max_apbs_points,
        ) if apbs is not None else ([], [])
        calculated_points = mask_point_cloud(calculated, h5f, args.max_label_points)
        provided_points = mask_point_cloud(provided, h5f, args.max_label_points)

        meta = (
            f"H5={h5_path.name} | case={case_name} | "
            f"APBS={apbs_name or 'none'} | label={provided_name or 'none'} | "
            f"box={h5f.attrs.get('box_size', '?')} | resolution={float(h5f.attrs.get('resolution', np.nan)):.3g} A/voxel"
        )

    warnings = []
    if protein_path is None or not protein_path.exists():
        warnings.append("Protein source file was not found; molecular protein model is unavailable.")
    if ligand_path is None or not ligand_path.exists():
        warnings.append("Ligand source file was not found; ligand sticks are unavailable.")
    warning_html = "".join(f"<div class='warn'>{warning}</div>" for warning in warnings)

    payload = {
        "title": f"{case_name} molecular feature viewer",
        "meta": meta,
        "proteinData": read_text(protein_path),
        "proteinFormat": structure_format(protein_path),
        "ligandData": read_text(ligand_path),
        "ligandFormat": structure_format(ligand_path),
        "apbsPositive": apbs_positive,
        "apbsNegative": apbs_negative,
        "calculatedLabel": calculated_points,
        "providedLabel": provided_points,
        "sphereRadius": args.sphere_radius,
        "warning_html": warning_html,
    }

    output = Path(args.output).expanduser().resolve() if args.output else h5_path.with_name(f"{h5_path.stem}_molecular_viewer.html")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(build_html(payload), encoding="utf-8")
    print(f"Protein: {protein_path if protein_path else 'not found'}")
    print(f"Ligand: {ligand_path if ligand_path else 'not found'}")
    print(f"Viewer: {output}")


if __name__ == "__main__":
    main()
