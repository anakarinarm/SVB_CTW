#!/usr/bin/env python3
"""Check MITgcm-written partial cells for the smaller-SVB smoke test."""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from MITgcmutils import mds
from scipy.ndimage import label
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


ROOT = Path(__file__).resolve().parents[1] / "ctw_reruns_telescoped_west64_south72"
RUN = Path("/Volumes/LaCie/CTW_SVB_reruns_telescoped_west64_south72/smoke_Rplus_50L_hfac")
OUT = ROOT / "smoke_Rplus_50L_hfac" / "hfac_review"


def main():
    OUT.mkdir(exist_ok=True)
    lon, lat, depth = [mds.rdmds(str(RUN / name)) for name in ("XC", "YC", "Depth")]
    dz = mds.rdmds(str(RUN / "DRF")).ravel()
    z = np.abs(mds.rdmds(str(RUN / "RC")).ravel())
    old = np.fromfile(ROOT / "west_south_input/bathy_realistic_plusSVB_W64S72.bin", dtype=">f8").reshape(depth.shape)
    new = np.fromfile(ROOT / "west_south_input/bathy_realistic_SVB50L_mouth_W64S72.bin", dtype=">f8").reshape(depth.shape)
    changed = old != new
    rows, cols = np.nonzero(changed)
    if not len(rows):
        raise RuntimeError("No bathymetry differences from Rplus")
    ys = slice(max(0, int(rows.min()) - 20), min(depth.shape[0], int(rows.max()) + 21))
    xs = slice(max(0, int(cols.min()) - 20), min(depth.shape[1], int(cols.max()) + 21))
    k = int(np.argmin(np.abs(z - 500)))
    fields = {}
    surface_faces = {}
    stats = {}
    for name in ("hFacC", "hFacW", "hFacS"):
        a = mds.rdmds(str(RUN / name))
        assert a.shape == (100, *depth.shape), (name, a.shape)
        assert np.all(np.isfinite(a)) and np.all((a >= 0) & (a <= 1)), name
        minimum = np.maximum(.1, np.minimum(2 / dz, 1))[:, None, None]
        assert not np.any((a > 0) & (a < minimum - 1e-12)), name
        stats[name] = dict(shape=list(a.shape), minimum=float(a.min()), maximum=float(a.max()),
                           partial_cells=int(np.count_nonzero((a > 0) & (a < 1))),
                           smallest_wet_fraction=float(a[a > 0].min()))
        fields[name] = (a[0, ys, xs].copy(), a[k, ys, xs].copy())
        surface_faces[name] = a[0].copy()
        if name == "hFacC":
            reconstructed = np.einsum("k,kji->ji", dz, a)
            assert np.allclose(reconstructed, depth, rtol=0, atol=1e-8)
            wet = a[0] > 0
            components, count = label(wet)
            sizes = np.bincount(components.ravel())
            ocean = int(sizes[1:].argmax() + 1)
            bay_wet = wet & changed
            isolated = wet & (components != ocean)
            stats["connectivity"] = dict(surface_components=int(count),
                                         wet_cells_outside_main_ocean=int(isolated.sum()),
                                         changed_region_wet_cells=int(bay_wet.sum()),
                                         changed_region_wet_cells_disconnected=int(np.count_nonzero(bay_wet & isolated)))
            stats["depth_reconstruction_max_error_m"] = float(np.max(np.abs(reconstructed - depth)))
        del a
    # Check actual open velocity faces, not only neighboring wet centers.
    node = np.full(wet.shape, -1, dtype=int)
    node[wet] = np.arange(wet.sum())
    ew = wet[:, :-1] & wet[:, 1:] & (surface_faces["hFacW"][:, 1:] > 0)
    ns = wet[:-1, :] & wet[1:, :] & (surface_faces["hFacS"][1:, :] > 0)
    source = np.concatenate((node[:, :-1][ew], node[:-1, :][ns]))
    target = np.concatenate((node[:, 1:][ew], node[1:, :][ns]))
    graph = coo_matrix((np.ones(len(source)), (source, target)), shape=(int(wet.sum()), int(wet.sum()))).tocsr()
    count, groups = connected_components(graph, directed=False)
    main_ocean = int(np.bincount(groups).argmax())
    stats["surface_face_connectivity"] = dict(
        components=int(count),
        wet_cells_disconnected_by_closed_faces=int(np.count_nonzero(groups != main_ocean)),
        changed_region_wet_cells_disconnected=int(np.count_nonzero(groups[node[bay_wet]] != main_ocean)))
    stats["input_depth_vs_processed_depth_max_difference_m"] = float(np.max(np.abs(-new - depth)))
    stats["changed_bathymetry_cells"] = int(changed.sum())
    stats["cells_dried_relative_to_Rplus"] = int(np.count_nonzero((old < 0) & (new == 0)))
    stats["cells_wetted_relative_to_Rplus"] = int(np.count_nonzero((old == 0) & (new < 0)))
    stats["plot_crop"] = dict(x_start=xs.start, x_end=xs.stop, y_start=ys.start, y_end=ys.stop)
    stats["deep_plot_level_1based"] = k + 1
    stats["deep_plot_depth_m"] = float(z[k])
    (OUT / "hfac_checks.json").write_text(json.dumps(stats, indent=2) + "\n")
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True, sharex=True, sharey=True)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_under("white")
    for col, (name, slices) in enumerate(fields.items()):
        for row, values in enumerate(slices):
            ax = axes[row, col]
            image = ax.pcolormesh(lon[ys, xs], lat[ys, xs], values, shading="auto", cmap=cmap, vmin=1e-10, vmax=1)
            ax.contour(lon[ys, xs], lat[ys, xs], old[ys, xs] < 0, levels=[.5], colors=".5", linewidths=.6, linestyles="--")
            ax.set_title(f"{name}: {'surface' if row == 0 else f'{z[k]:.1f} m'}")
            ax.set(xlabel="Longitude (degrees E)", ylabel="Latitude (degrees N)", aspect="equal")
    fig.colorbar(image, ax=axes.ravel().tolist(), label="Wet fraction (white = dry)", shrink=.85)
    fig.suptitle("Rplus_50L partial cells; dashed outline = original Rplus coastline")
    fig.savefig(OUT / "hfac_bay_inspection.png", dpi=160)
    fig.savefig(OUT / "hfac_bay_inspection.pdf")
    plt.close(fig)
    print(json.dumps(stats, indent=2), flush=True)
    print(f"Wrote {OUT / 'hfac_bay_inspection.png'}", flush=True)


if __name__ == "__main__":
    main()
