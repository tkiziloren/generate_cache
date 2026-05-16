#!/usr/bin/env python3
"""Visualize H5 cache feature and label channels.

The script is intentionally file-oriented: it accepts one or more H5 files and
writes static PNG summaries that are easy to inspect locally or attach to notes.
It can also write an interactive Plotly HTML page containing every feature and
label channel, including future v2 feature names.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
except Exception:  # pragma: no cover - static PNG output does not need Plotly.
    make_subplots = None
    go = None


DENSE_CONTINUOUS_FEATURES = {
    "electrostatic_grid",
    "dist_to_surface",
    "dist_to_ligand",
    "hydrophobicity",
}

SIGNED_FEATURE_HINTS = {
    "electrostatic",
    "charge",
    "partialcharge",
    "hydrophobicity",
}

DEFAULT_SAMPLE_DIR = Path(
    "/Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/codon_h5_samples_2026-05-16"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Visualize all feature and label channels in H5 cache files.")
    parser.add_argument("h5_paths", nargs="*", help="One or more .h5 files.")
    parser.add_argument("--input-dir", default=None, help="Directory to search for H5 files.")
    parser.add_argument("--glob", default="*box36*.h5", help="Glob used with --input-dir.")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for generated figures. Defaults to <input parent>/visualizations.",
    )
    parser.add_argument("--cols", type=int, default=6, help="Number of columns in the PNG grid.")
    parser.add_argument("--dpi", type=int, default=150, help="PNG DPI.")
    parser.add_argument("--histograms", action="store_true", help="Also write histogram PNGs.")
    parser.add_argument("--html", action="store_true", help="Also write interactive Plotly 3D HTML pages.")
    parser.add_argument(
        "--presentation-html",
        action="store_true",
        help="Write the curated presentation-style 3D HTML layout used in earlier thesis figures.",
    )
    parser.add_argument("--html-cols", type=int, default=4, help="Number of 3D panels per HTML row.")
    parser.add_argument("--html-max-points", type=int, default=5000, help="Max scatter points per sparse 3D channel.")
    parser.add_argument(
        "--html-self-contained",
        action="store_true",
        help="Embed Plotly JS inside each HTML file instead of loading it from CDN.",
    )
    return parser.parse_args()


def discover_h5_paths(args: argparse.Namespace) -> list[Path]:
    paths = [Path(path).expanduser().resolve() for path in args.h5_paths]
    if args.input_dir:
        paths.extend(sorted(Path(args.input_dir).expanduser().resolve().rglob(args.glob)))
    if not paths and DEFAULT_SAMPLE_DIR.exists():
        paths = sorted(DEFAULT_SAMPLE_DIR.rglob(args.glob))

    unique_paths = []
    seen = set()
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        if not path.exists():
            raise FileNotFoundError(path)
        if path.suffix != ".h5":
            raise ValueError(f"Expected .h5 file: {path}")
        unique_paths.append(path)
    if not unique_paths:
        raise SystemExit("No H5 files found. Pass paths or use --input-dir.")
    return unique_paths


def output_dir_for(args: argparse.Namespace, first_h5: Path) -> Path:
    if args.output_dir:
        out_dir = Path(args.output_dir).expanduser().resolve()
    else:
        out_dir = first_h5.parent / "visualizations"
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def read_channel_names(h5f: h5py.File) -> list[tuple[str, str]]:
    channels = []
    if "features" in h5f:
        channels.extend(("features", name) for name in sorted(h5f["features"].keys()))
    if "auxiliary" in h5f:
        channels.extend(("auxiliary", name) for name in sorted(h5f["auxiliary"].keys()))
    if "label" in h5f:
        channels.extend(("label", name) for name in sorted(h5f["label"].keys()))
    return channels


def classify_channel(group: str, name: str) -> str:
    if group == "label":
        return "etiketler"
    if group == "auxiliary":
        return "yardimci"
    if name.startswith("atomic_"):
        return "atomic"
    if name.startswith("electrostatic_grid_v1") or name == "electrostatic_grid":
        return "apbs_v1"
    if name.startswith("electrostatic_grid_v2"):
        return "apbs_v2"
    if name in {"shape", "hydrophobicity", "dist_to_surface", "dist_to_ligand", "ligand"}:
        return "yardimci"
    return "diger"


def grouped_channels(channels: list[tuple[str, str]]) -> list[tuple[str, list[tuple[str, str]]]]:
    order = ["atomic", "apbs_v1", "apbs_v2", "yardimci", "etiketler", "diger", "tum_kanallar"]
    always_visible = {"atomic", "apbs_v1", "apbs_v2", "yardimci"}
    groups = {name: [] for name in order}
    for group, name in channels:
        groups.setdefault(classify_channel(group, name), []).append((group, name))
    groups["tum_kanallar"] = channels
    return [(name, groups[name]) for name in order if groups.get(name) or name in always_visible]


def read_named_array(h5f: h5py.File, name: str, groups: tuple[str, ...] = ("features", "auxiliary", "label")):
    for group in groups:
        path = f"{group}/{name}"
        if path in h5f:
            return h5f[path][:]
    return None


def first_available_name(h5f: h5py.File, names: list[str], groups: tuple[str, ...] = ("features", "auxiliary", "label")):
    for name in names:
        for group in groups:
            if f"{group}/{name}" in h5f:
                return group, name
    return None, None


def is_binary_like(data: np.ndarray) -> bool:
    if data.dtype == np.bool_:
        return True
    if np.issubdtype(data.dtype, np.integer):
        finite_values = np.unique(data)
        return finite_values.size <= 2 and set(finite_values.tolist()).issubset({0, 1})
    if data.size > 250_000:
        sample = data.ravel()[:: max(1, data.size // 250_000)]
    else:
        sample = data.ravel()
    finite_values = np.unique(sample[np.isfinite(sample)])
    return finite_values.size <= 2 and set(np.round(finite_values, 6).tolist()).issubset({0.0, 1.0})


def nonzero_fraction(data: np.ndarray) -> float:
    if data.size == 0:
        return 0.0
    return float(np.count_nonzero(data)) / float(data.size)


def choose_focus_slice(h5f: h5py.File) -> int:
    candidates = []
    if "label" in h5f:
        for name in h5f["label"]:
            data = h5f[f"label/{name}"][:]
            if data.ndim == 3:
                candidates.append(data > 0)
    if "features/ligand" in h5f:
        candidates.append(h5f["features/ligand"][:] > 0)
    if "auxiliary/ligand" in h5f:
        candidates.append(h5f["auxiliary/ligand"][:] > 0)
    if "features/shape" in h5f:
        candidates.append(h5f["features/shape"][:] > 0)
    if not candidates:
        first_channel = h5f[read_channel_names(h5f)[0][0]][read_channel_names(h5f)[0][1]]
        return int(first_channel.shape[-1] // 2)

    summed = np.zeros_like(candidates[0], dtype=np.int32)
    for mask in candidates:
        summed += mask.astype(np.int32)
    z_profile = summed.sum(axis=(0, 1))
    if z_profile.max() == 0:
        return int(summed.shape[-1] // 2)
    return int(np.argmax(z_profile))


def should_use_projection(group: str, name: str, data: np.ndarray) -> bool:
    if group == "label":
        return True
    if is_binary_like(data):
        return True
    if name.startswith("atomic_"):
        return True
    if nonzero_fraction(data) < 0.05:
        return True
    return False


def panel_image(group: str, name: str, data: np.ndarray, focus_z: int) -> tuple[np.ndarray, str]:
    if data.ndim != 3:
        raise ValueError(f"Expected 3D channel for {group}/{name}, got {data.shape}")
    if should_use_projection(group, name, data):
        if np.nanmin(data) < 0 and np.nanmax(data) > 0:
            idx = np.nanargmax(np.abs(data), axis=2)
            image = np.take_along_axis(data, idx[:, :, None], axis=2)[:, :, 0]
            mode = "maxabs-z"
        else:
            image = np.nanmax(data, axis=2)
            mode = "max-z"
    else:
        focus_z = min(max(0, focus_z), data.shape[2] - 1)
        image = data[:, :, focus_z]
        mode = f"z={focus_z}"
    return np.rot90(image), mode


def color_limits(name: str, data: np.ndarray) -> tuple[float | None, float | None, str]:
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return None, None, "viridis"
    lname = name.lower()
    signed = any(hint in lname for hint in SIGNED_FEATURE_HINTS) and finite.min() < 0 < finite.max()
    if signed:
        vmax = float(np.percentile(np.abs(finite), 99.0))
        vmax = vmax if vmax > 0 else 1.0
        return -vmax, vmax, "RdBu_r"
    if is_binary_like(data):
        return 0.0, 1.0, "gray_r"
    vmin = float(np.percentile(finite, 1.0))
    vmax = float(np.percentile(finite, 99.0))
    if math.isclose(vmin, vmax):
        vmin = float(finite.min())
        vmax = float(finite.max())
    if math.isclose(vmin, vmax):
        return None, None, "viridis"
    return vmin, vmax, "viridis"


def channel_stats(data: np.ndarray) -> str:
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return "no finite values"
    nonzero = int(np.count_nonzero(data))
    return (
        f"min {finite.min():.2g} | max {finite.max():.2g} | "
        f"sum {np.nansum(data):.2g} | nz {nonzero}"
    )


def write_channel_grid_png(h5_path: Path, out_dir: Path, cols: int, dpi: int) -> Path:
    with h5py.File(h5_path, "r") as h5f:
        channels = read_channel_names(h5f)
        focus_z = choose_focus_slice(h5f)
        if not channels:
            raise RuntimeError(f"No feature/label channels found in {h5_path}")

        rows = math.ceil(len(channels) / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(cols * 3.2, rows * 3.25), squeeze=False)
        fig.suptitle(
            (
                f"{h5_path.name} | box={h5f.attrs.get('box_size', '?')} | "
                f"resolution={float(h5f.attrs.get('resolution', np.nan)):.3g} A/voxel | "
                f"span={float(h5f.attrs.get('physical_span_angstrom', np.nan)):.3g} A"
            ),
            fontsize=14,
            fontweight="bold",
        )

        for ax in axes.ravel():
            ax.axis("off")

        for ax, (group, name) in zip(axes.ravel(), channels):
            data = h5f[f"{group}/{name}"][:]
            image, mode = panel_image(group, name, data, focus_z)
            vmin, vmax, cmap = color_limits(name, data)
            ax.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
            ax.set_title(f"{group}/{name}\n{mode} | {channel_stats(data)}", fontsize=7)
            ax.axis("off")

        fig.tight_layout(rect=(0, 0, 1, 0.965))
        out_path = out_dir / f"{h5_path.stem}_all_channels.png"
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return out_path


def write_histogram_png(h5_path: Path, out_dir: Path, cols: int, dpi: int) -> Path:
    with h5py.File(h5_path, "r") as h5f:
        channels = read_channel_names(h5f)
        rows = math.ceil(len(channels) / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(cols * 3.0, rows * 2.4), squeeze=False)
        fig.suptitle(f"{h5_path.name} channel histograms", fontsize=14, fontweight="bold")
        for ax in axes.ravel():
            ax.axis("off")
        for ax, (group, name) in zip(axes.ravel(), channels):
            data = h5f[f"{group}/{name}"][:].ravel()
            finite = data[np.isfinite(data)]
            if finite.size == 0:
                continue
            if finite.size > 250_000:
                finite = finite[:: max(1, finite.size // 250_000)]
            ax.hist(finite, bins=50, color="#3b82f6", alpha=0.75)
            ax.set_title(f"{group}/{name}", fontsize=7)
            ax.tick_params(labelsize=6)
            ax.axis("on")
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        out_path = out_dir / f"{h5_path.stem}_histograms.png"
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return out_path


def sparse_points(data: np.ndarray, max_points: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    coords = np.argwhere(data > 0)
    if coords.size == 0:
        return np.array([]), np.array([]), np.array([])
    if len(coords) > max_points:
        step = math.ceil(len(coords) / max_points)
        coords = coords[::step]
    return coords[:, 0], coords[:, 1], coords[:, 2]


def sparse_value_points(data: np.ndarray, max_points: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if is_binary_like(data):
        coords = np.argwhere(data > 0)
    else:
        coords = np.argwhere(np.isfinite(data) & (np.abs(data) > 1.0e-8))
    if coords.size == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])
    if len(coords) > max_points:
        step = math.ceil(len(coords) / max_points)
        coords = coords[::step]
    values = data[coords[:, 0], coords[:, 1], coords[:, 2]]
    return coords[:, 0], coords[:, 1], coords[:, 2], values


def grid_coordinates(shape: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.arange(shape[0]).repeat(shape[1] * shape[2])
    y = np.tile(np.arange(shape[1]).repeat(shape[2]), shape[0])
    z = np.tile(np.arange(shape[2]), shape[0] * shape[1])
    return x, y, z


def channel_kind(group: str, name: str, data: np.ndarray) -> str:
    if group == "label" or name in {"shape", "ligand"} or is_binary_like(data):
        return "mask"
    if name.startswith("atomic_") or nonzero_fraction(data) < 0.08:
        return "sparse"
    return "dense"


def sparse_color_scale(name: str, values: np.ndarray):
    if values.size == 0:
        return "#7c3aed", None, None
    signed = np.nanmin(values) < 0 < np.nanmax(values)
    if signed or any(hint in name.lower() for hint in SIGNED_FEATURE_HINTS):
        vmax = float(np.nanpercentile(np.abs(values), 99.0))
        vmax = vmax if vmax > 0 else 1.0
        return values, -vmax, vmax
    return values, None, None


def add_channel_trace(fig, data: np.ndarray, group: str, name: str, row: int, col: int, max_points: int):
    full_name = f"{group}/{name}"
    kind = channel_kind(group, name, data)
    box_size = data.shape[0]

    if kind == "dense":
        finite = data[np.isfinite(data)]
        if finite.size == 0 or math.isclose(float(finite.min()), float(finite.max())):
            fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode="text", text=[f"{full_name}: constant"]), row=row, col=col)
            return
        x, y, z = grid_coordinates(data.shape)
        lname = name.lower()
        signed = any(hint in lname for hint in SIGNED_FEATURE_HINTS) and finite.min() < 0 < finite.max()
        if signed:
            vmax = float(np.percentile(np.abs(finite), 97.5))
            vmax = vmax if vmax > 0 else float(np.max(np.abs(finite)))
            isomin, isomax = -vmax, vmax
            colorscale = "RdBu_r"
        else:
            isomin = float(np.percentile(finite, 70.0))
            isomax = float(np.percentile(finite, 99.0))
            if math.isclose(isomin, isomax):
                isomin, isomax = float(finite.min()), float(finite.max())
            colorscale = "Viridis"
        fig.add_trace(
            go.Isosurface(
                x=x,
                y=y,
                z=z,
                value=data.ravel(),
                isomin=isomin,
                isomax=isomax,
                surface_count=3,
                opacity=0.55,
                colorscale=colorscale,
                caps=dict(x_show=False, y_show=False, z_show=False),
                showscale=False,
                name=full_name,
                showlegend=False,
            ),
            row=row,
            col=col,
        )
        return

    if kind == "mask":
        x, y, z = sparse_points(data, max_points=max_points)
        if x.size == 0:
            fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode="text", text=[f"{full_name}: empty"]), row=row, col=col)
            return
        color = "#2563eb"
        if group == "label":
            color = "#dc2626"
        elif name == "ligand":
            color = "#16a34a"
        elif name == "shape":
            color = "#64748b"
        fig.add_trace(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="markers",
                marker=dict(size=2.4, color=color, opacity=0.72),
                name=full_name,
                showlegend=False,
            ),
            row=row,
            col=col,
        )
        return

    x, y, z, values = sparse_value_points(data, max_points=max_points)
    if x.size == 0:
        fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode="text", text=[f"{full_name}: empty"]), row=row, col=col)
        return
    marker_color, cmin, cmax = sparse_color_scale(name, values)
    marker = dict(size=2.4, color=marker_color, opacity=0.78, colorscale="RdBu_r" if cmin is not None else "Viridis")
    if cmin is not None:
        marker["cmin"] = cmin
        marker["cmax"] = cmax
    fig.add_trace(
        go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="markers",
            marker=marker,
            name=full_name,
            showlegend=False,
        ),
        row=row,
        col=col,
    )


def add_isosurface_trace(
    fig,
    data: np.ndarray,
    row: int,
    col: int,
    name: str,
    isomin: float,
    isomax: float,
    colorscale: str,
    opacity: float,
    surface_count: int = 1,
    showscale: bool = False,
) -> None:
    x, y, z = grid_coordinates(data.shape)
    fig.add_trace(
        go.Isosurface(
            x=x,
            y=y,
            z=z,
            value=data.ravel(),
            isomin=isomin,
            isomax=isomax,
            surface_count=surface_count,
            opacity=opacity,
            colorscale=colorscale,
            caps=dict(x_show=False, y_show=False, z_show=False),
            showscale=showscale,
            name=name,
            showlegend=False,
        ),
        row=row,
        col=col,
    )


def add_mask_points(
    fig,
    data: np.ndarray,
    row: int,
    col: int,
    name: str,
    color: str,
    max_points: int,
    size: float = 2.6,
    opacity: float = 0.75,
) -> None:
    x, y, z = sparse_points(data, max_points=max_points)
    if x.size == 0:
        fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode="text", text=["No Data"]), row=row, col=col)
        return
    fig.add_trace(
        go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="markers",
            marker=dict(size=size, color=color, opacity=opacity),
            name=name,
            showlegend=False,
        ),
        row=row,
        col=col,
    )


def add_value_points(
    fig,
    data: np.ndarray,
    row: int,
    col: int,
    name: str,
    max_points: int,
    size: float = 2.6,
    opacity: float = 0.78,
) -> None:
    x, y, z, values = sparse_value_points(data, max_points=max_points)
    if x.size == 0:
        fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode="text", text=["No Data"]), row=row, col=col)
        return
    color, cmin, cmax = sparse_color_scale(name, values)
    marker = dict(size=size, color=color, opacity=opacity, colorscale="RdBu_r" if cmin is not None else "Viridis")
    if cmin is not None:
        marker["cmin"] = cmin
        marker["cmax"] = cmax
    fig.add_trace(
        go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="markers",
            marker=marker,
            name=name,
            showlegend=False,
        ),
        row=row,
        col=col,
    )


def add_surface_binding_panel(
    fig,
    shape: np.ndarray,
    label: np.ndarray | None,
    ligand: np.ndarray | None,
    row: int,
    col: int,
    max_points: int,
) -> None:
    add_isosurface_trace(
        fig,
        shape.astype(np.float32),
        row,
        col,
        "protein_surface",
        isomin=0.5,
        isomax=1.0,
        colorscale="Blues",
        opacity=0.15,
    )
    if label is not None:
        add_mask_points(fig, label, row, col, "binding_site", "#dc2626", max_points=max_points, size=2.8, opacity=0.9)
    if ligand is not None:
        add_mask_points(fig, ligand, row, col, "ligand", "#16a34a", max_points=max_points, size=3.2, opacity=0.9)


def choose_presentation_apbs(h5f: h5py.File) -> tuple[str | None, np.ndarray | None]:
    candidates = [
        "electrostatic_grid_v1_ligand_proximal_chains_7A_raw",
        "electrostatic_grid",
        "electrostatic_grid_v2_full_protein_raw",
    ]
    group, name = first_available_name(h5f, candidates, groups=("features", "auxiliary"))
    if name is None:
        return None, None
    return name, h5f[f"{group}/{name}"][:]


def choose_provided_label(h5f: h5py.File) -> tuple[str | None, np.ndarray | None]:
    candidates = [
        "binding_site_in_dataset",
        "binding_site_fpocket_selected",
        "binding_site_cavity6",
        "binding_site_cavityALL",
        "binding_site_site",
    ]
    group, name = first_available_name(h5f, candidates, groups=("label",))
    if name is None:
        return None, None
    return name, h5f[f"{group}/{name}"][:]


def write_presentation_html(
    h5_path: Path,
    out_dir: Path,
    max_points: int,
    self_contained: bool,
    apbs_min: float = -15.0,
    apbs_max: float = 15.0,
) -> Path | None:
    if make_subplots is None or go is None:
        print("[WARN] Plotly is not available; skipping presentation HTML output.")
        return None

    with h5py.File(h5_path, "r") as h5f:
        atomic_features = sorted(name for name in h5f.get("features", {}) if name.startswith("atomic_"))
        key_titles = [
            "APBS",
            "Yüzeye uzaklık",
            "Liganda uzaklık",
            "Protein şekli",
            "Hesaplanan cep ile protein",
            "Veri seti cebi ile protein",
        ]
        ncols = 6
        nrows = 1 + math.ceil(len(atomic_features) / ncols)
        subplot_titles = key_titles + atomic_features
        fig = make_subplots(
            rows=nrows,
            cols=ncols,
            specs=[[{"type": "scene"} for _ in range(ncols)] for _ in range(nrows)],
            subplot_titles=subplot_titles,
            horizontal_spacing=0.025,
            vertical_spacing=0.035,
        )

        apbs_name, apbs = choose_presentation_apbs(h5f)
        shape = read_named_array(h5f, "shape", groups=("features",))
        ligand = read_named_array(h5f, "ligand", groups=("features", "auxiliary"))
        dist_to_surface = read_named_array(h5f, "dist_to_surface", groups=("features", "auxiliary"))
        dist_to_ligand = read_named_array(h5f, "dist_to_ligand", groups=("features", "auxiliary"))
        calculated = read_named_array(h5f, "binding_site_calculated", groups=("label",))
        provided_name, provided = choose_provided_label(h5f)

        if apbs is not None:
            add_isosurface_trace(
                fig,
                apbs.astype(np.float32),
                1,
                1,
                apbs_name or "APBS",
                isomin=apbs_min,
                isomax=apbs_max,
                colorscale="RdBu",
                opacity=0.55,
                surface_count=2,
                showscale=True,
            )
        if dist_to_surface is not None:
            finite = dist_to_surface[np.isfinite(dist_to_surface)]
            nonzero = finite[np.abs(finite) > 1.0e-8]
            iso = float(np.percentile(nonzero, 98)) if nonzero.size else 0.1
            add_isosurface_trace(
                fig,
                dist_to_surface.astype(np.float32),
                1,
                2,
                "dist_to_surface",
                isomin=iso,
                isomax=float(np.max(finite)) if finite.size else 1.0,
                colorscale="Viridis",
                opacity=0.70,
                surface_count=2,
            )
        if dist_to_ligand is not None:
            finite = dist_to_ligand[np.isfinite(dist_to_ligand)]
            nonzero = finite[np.abs(finite) > 1.0e-8]
            iso = float(np.percentile(nonzero, 98)) if nonzero.size else 0.1
            add_isosurface_trace(
                fig,
                dist_to_ligand.astype(np.float32),
                1,
                3,
                "dist_to_ligand",
                isomin=iso,
                isomax=float(np.max(finite)) if finite.size else 1.0,
                colorscale="Magma",
                opacity=0.70,
                surface_count=2,
            )
        if shape is not None:
            add_isosurface_trace(
                fig,
                shape.astype(np.float32),
                1,
                4,
                "shape",
                isomin=0.5,
                isomax=1.0,
                colorscale="Blues",
                opacity=0.18,
            )
            add_surface_binding_panel(fig, shape, calculated, ligand, 1, 5, max_points=max_points)
            add_surface_binding_panel(fig, shape, provided, ligand, 1, 6, max_points=max_points)

        for idx, feature_name in enumerate(atomic_features):
            row = idx // ncols + 2
            col = idx % ncols + 1
            data = h5f[f"features/{feature_name}"][:]
            add_value_points(fig, data, row, col, feature_name, max_points=max_points)

        box_size = int(h5f.attrs.get("box_size", shape.shape[0] if shape is not None else 36))
        for idx in range(1, len(subplot_titles) + 1):
            scene_name = f"scene{idx}" if idx > 1 else "scene"
            fig.update_layout(
                {
                    scene_name: dict(
                        aspectmode="cube",
                        xaxis=dict(range=[0, box_size], title="x"),
                        yaxis=dict(range=[0, box_size], title="y"),
                        zaxis=dict(range=[0, box_size], title="z"),
                    )
                }
            )

        fig.update_layout(
            title=dict(
                text=(
                    "Yapılan Diğer Çalışmalar - Özniteliklerin Görselleştirilmesi"
                    f"<br><sup>{h5_path.name}"
                    f" | box={h5f.attrs.get('box_size', '?')}"
                    f" | çözünürlük={float(h5f.attrs.get('resolution', np.nan)):.3g} A/voxel"
                    f" | APBS={apbs_name or 'yok'}"
                    f" | dataset label={provided_name or 'yok'}</sup>"
                ),
                x=0.01,
                xanchor="left",
                font=dict(size=30, color="#111827"),
            ),
            width=2100,
            height=max(980, 310 * nrows + 190),
            showlegend=False,
            paper_bgcolor="#f2f2f2",
            plot_bgcolor="#ffffff",
            margin=dict(l=35, r=35, t=135, b=35),
        )
        fig.update_annotations(font_size=11)
        out_path = out_dir / f"{h5_path.stem}_presentation_3d.html"
        fig.write_html(out_path, include_plotlyjs=True if self_contained else "cdn")
        return out_path


def write_html_overview(
    h5_path: Path,
    out_dir: Path,
    max_points: int,
    cols: int,
    self_contained: bool,
) -> Path | None:
    if make_subplots is None or go is None:
        print("[WARN] Plotly is not available; skipping HTML output.")
        return None

    with h5py.File(h5_path, "r") as h5f:
        channels = read_channel_names(h5f)
        if not channels:
            return None
        first_group, first_name = channels[0]
        box_size = int(h5f.attrs.get("box_size", h5f[first_group][first_name].shape[0]))

        tab_parts = []
        buttons = []
        include_plotlyjs = True if self_contained else "cdn"
        plotly_js_written = False
        for tab_index, (tab_name, tab_channels) in enumerate(grouped_channels(channels)):
            if tab_channels:
                tab_cols = max(1, min(cols, len(tab_channels)))
                rows = math.ceil(len(tab_channels) / tab_cols)
                fig = make_subplots(
                    rows=rows,
                    cols=tab_cols,
                    specs=[[{"type": "scene"} for _ in range(tab_cols)] for _ in range(rows)],
                    subplot_titles=[f"{group}/{name}" for group, name in tab_channels],
                    horizontal_spacing=0.015,
                    vertical_spacing=0.055,
                )
                for idx, (group, name) in enumerate(tab_channels):
                    row = idx // tab_cols + 1
                    col = idx % tab_cols + 1
                    data = h5f[f"{group}/{name}"][:]
                    add_channel_trace(fig, data, group, name, row, col, max_points)

                for idx in range(1, len(tab_channels) + 1):
                    scene_name = f"scene{idx}" if idx > 1 else "scene"
                    fig.update_layout(
                        {
                            scene_name: dict(
                                aspectmode="cube",
                                xaxis=dict(range=[0, box_size], title="x"),
                                yaxis=dict(range=[0, box_size], title="y"),
                                zaxis=dict(range=[0, box_size], title="z"),
                            )
                        }
                    )
                fig.update_layout(
                    title=(
                        f"{h5_path.name} | {tab_name} | "
                        f"box={h5f.attrs.get('box_size', '?')} | "
                        f"resolution={float(h5f.attrs.get('resolution', np.nan)):.3g} A/voxel"
                    ),
                    width=max(1000, 380 * tab_cols),
                    height=max(650, 390 * rows),
                    showlegend=False,
                    margin=dict(l=10, r=10, t=72, b=10),
                )
                fig.update_annotations(font_size=10)
                div = fig.to_html(
                    full_html=False,
                    include_plotlyjs=include_plotlyjs if not plotly_js_written else False,
                )
                plotly_js_written = True
            else:
                div = (
                    "<div class='empty-tab'>"
                    f"<h2>{tab_name}</h2>"
                    "<p>Bu H5 dosyasında bu gruba ait kanal yok.</p>"
                    "</div>"
                )
            display = "block" if tab_index == 0 else "none"
            tab_id = f"tab-{tab_name}"
            tab_parts.append(f"<section id='{tab_id}' class='tab-panel' style='display:{display}'>{div}</section>")
            active = " active" if tab_index == 0 else ""
            buttons.append(
                f"<button class='tab-button{active}' data-target='{tab_id}'>{tab_name}"
                f"<span>{len(tab_channels)}</span></button>"
            )

        out_path = out_dir / f"{h5_path.stem}_all_channels_3d.html"
        html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{h5_path.name} feature viewer</title>
  <style>
    body {{ margin:0; font-family:Arial, sans-serif; background:#f8fafc; color:#0f172a; }}
    header {{ padding:18px 24px 10px; background:#ffffff; border-bottom:1px solid #e2e8f0; position:sticky; top:0; z-index:20; }}
    h1 {{ margin:0 0 10px; font-size:22px; }}
    .meta {{ color:#475569; font-size:13px; }}
    .tabs {{ display:flex; gap:8px; flex-wrap:wrap; margin-top:14px; }}
    .tab-button {{ border:1px solid #cbd5e1; background:#fff; color:#0f172a; padding:8px 12px; border-radius:6px; cursor:pointer; font-size:13px; }}
    .tab-button.active {{ background:#0f172a; color:#fff; border-color:#0f172a; }}
    .tab-button span {{ margin-left:8px; opacity:.75; font-size:11px; }}
    main {{ padding:16px; }}
    .tab-panel {{ background:#fff; border:1px solid #e2e8f0; border-radius:8px; overflow:auto; }}
    .empty-tab {{ min-height:260px; display:flex; flex-direction:column; align-items:center; justify-content:center; color:#475569; }}
    .empty-tab h2 {{ margin:0 0 8px; color:#0f172a; font-size:22px; }}
    .empty-tab p {{ margin:0; font-size:14px; }}
  </style>
</head>
<body>
  <header>
    <h1>H5 Feature Viewer</h1>
    <div class="meta">{h5_path.name} | box={h5f.attrs.get('box_size', '?')} | resolution={float(h5f.attrs.get('resolution', np.nan)):.3g} A/voxel</div>
    <nav class="tabs">{''.join(buttons)}</nav>
  </header>
  <main>{''.join(tab_parts)}</main>
  <script>
    document.querySelectorAll('.tab-button').forEach((button) => {{
      button.addEventListener('click', () => {{
        document.querySelectorAll('.tab-button').forEach((b) => b.classList.remove('active'));
        document.querySelectorAll('.tab-panel').forEach((panel) => panel.style.display = 'none');
        button.classList.add('active');
        document.getElementById(button.dataset.target).style.display = 'block';
        window.dispatchEvent(new Event('resize'));
      }});
    }});
  </script>
</body>
</html>
"""
        out_path.write_text(html, encoding="utf-8")
        return out_path


def write_index_html(out_dir: Path, html_paths: list[Path], png_paths: list[Path], presentation_paths: list[Path]) -> Path:
    index_path = out_dir / "index.html"
    lines = [
        "<!doctype html>",
        "<html><head><meta charset='utf-8'><title>H5 Feature Visualizations</title>",
        "<style>body{font-family:Arial,sans-serif;margin:24px;line-height:1.45}"
        "h1{font-size:24px} h2{font-size:18px;margin-top:24px}"
        "li{margin:6px 0} code{background:#f1f5f9;padding:2px 4px;border-radius:4px}</style>",
        "</head><body>",
        "<h1>H5 Feature Visualizations</h1>",
        "<p>Interactive pages support rotate, pan, zoom, and browser zoom. "
        "Every <code>features/*</code> and <code>label/*</code> channel present in the H5 file is included, "
        "including v2 feature names.</p>",
        "<h2>Presentation-Style 3D HTML</h2>",
        "<ul>",
    ]
    for path in presentation_paths:
        lines.append(f"<li><a href='{path.name}'>{path.name}</a></li>")
    lines.extend([
        "</ul>",
        "<h2>Interactive 3D HTML</h2>",
        "<ul>",
    ])
    for path in html_paths:
        lines.append(f"<li><a href='{path.name}'>{path.name}</a></li>")
    lines.extend(["</ul>", "<h2>Static PNG Summaries</h2>", "<ul>"])
    for path in png_paths:
        lines.append(f"<li><a href='{path.name}'>{path.name}</a></li>")
    lines.extend(["</ul>", "</body></html>"])
    index_path.write_text("\n".join(lines), encoding="utf-8")
    return index_path


def main() -> None:
    args = parse_args()
    h5_paths = discover_h5_paths(args)
    out_dir = output_dir_for(args, h5_paths[0])

    print(f"Output directory: {out_dir}")
    html_paths = []
    presentation_paths = []
    png_paths = []
    for h5_path in h5_paths:
        print(f"\n[H5] {h5_path}")
        grid_png = write_channel_grid_png(h5_path, out_dir, cols=args.cols, dpi=args.dpi)
        png_paths.append(grid_png)
        print(f"  wrote {grid_png}")
        if args.histograms:
            hist_png = write_histogram_png(h5_path, out_dir, cols=args.cols, dpi=args.dpi)
            png_paths.append(hist_png)
            print(f"  wrote {hist_png}")
        if args.html:
            html_path = write_html_overview(
                h5_path,
                out_dir,
                max_points=args.html_max_points,
                cols=args.html_cols,
                self_contained=args.html_self_contained,
            )
            if html_path:
                html_paths.append(html_path)
                print(f"  wrote {html_path}")
        if args.presentation_html:
            presentation_path = write_presentation_html(
                h5_path,
                out_dir,
                max_points=args.html_max_points,
                self_contained=args.html_self_contained,
            )
            if presentation_path:
                presentation_paths.append(presentation_path)
                print(f"  wrote {presentation_path}")
    if html_paths or png_paths or presentation_paths:
        index_path = write_index_html(
            out_dir,
            html_paths=html_paths,
            png_paths=png_paths,
            presentation_paths=presentation_paths,
        )
        print(f"\nIndex: {index_path}")


if __name__ == "__main__":
    main()
