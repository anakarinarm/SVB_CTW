#!/usr/bin/env python3
"""Create isolated 30-day realistic-bathymetry Gaussian-wind experiments."""
from __future__ import annotations

import argparse
import csv
import difflib
import json
import math
import re
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "ctw_reruns_telescoped_west64_south72"
INPUT = ROOT / "west_south_input"
NEW_INPUT = ROOT / "gaussian5d_input"
REPORT = ROOT / "gaussian5d_review"
EXTERNAL = Path("/Volumes/LaCie/CTW_SVB_reruns_telescoped_west64_south72")
NX, NY, DT, SPACING = 576, 684, 60, 1800
END, CYCLE = 30 * 86400, 31 * 86400
NREC = CYCLE // SPACING
TC, SIGMA, TAU = 216000., 72000., 0.1
FILES = ("data", "data.obcs", "data.pkg", "data.diagnostics", "eedata")


def gaussian(t):
    return np.exp(-0.5 * ((np.asarray(t) - TC) / SIGMA) ** 2)


def wind_name(component):
    return f"upwelling_{component}_G5d_31day_576x684_W64S72.bin"


def old_wind(component):
    p = INPUT / f"upwelling_{component}_halfgauss2hr_16day_576x684_W64S72.bin"
    return np.memmap(p, dtype=">f8", mode="r", shape=(768, NY, NX))


def interpolate(records, t, spacing=SPACING):
    q = np.asarray(t) / spacing - 0.5
    lower = np.floor(q).astype(int)
    weight = q - lower
    return records[lower % len(records)] * (1 - weight) + records[(lower + 1) % len(records)] * weight


def replace(text, name, value):
    pattern = rf"(?m)^(\s*{name}\s*=)[^,\n]*,"
    result, count = re.subn(pattern, lambda m: m[1] + value + ",", text)
    if count != 1:
        raise ValueError(f"Expected one assignment for {name}, found {count}")
    return result


def link_new(path, target):
    if path.is_symlink() and path.resolve() == target.resolve():
        return
    if path.exists() or path.is_symlink():
        raise FileExistsError(path)
    path.symlink_to(target.resolve())


def stage_external():
    if not EXTERNAL.parent.is_mount():
        raise RuntimeError("LaCie must be mounted")
    for base in ("Rplus", "Rminus"):
        out = EXTERNAL / (base + "_G5d")
        if out.exists():
            raise FileExistsError(f"Refusing to modify existing external run: {out}")
    for base in ("Rplus", "Rminus"):
        out = EXTERNAL / (base + "_G5d")
        out.mkdir()
        for name in FILES:
            link_new(out / name, ROOT / out.name / name)
        # Per-run links do not modify any existing family-level links.
        link_new(out / "mitgcmuv", ROOT / "build" / "mitgcmuv")
    link_new(EXTERNAL / "gaussian5d_input", NEW_INPUT)
    print(f"Staged new runs at {EXTERNAL}; no simulation launched", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-external", action="store_true")
    parser.add_argument("--review-only", action="store_true", help="Refresh comparison artifacts without changing inputs")
    args = parser.parse_args()
    if args.stage_external:
        stage_external()
        return
    if not args.review_only:
        for p in (NEW_INPUT, REPORT, ROOT / "Rplus_G5d", ROOT / "Rminus_G5d"):
            if p.exists():
                raise FileExistsError(f"Refusing to overwrite {p}")
        NEW_INPUT.mkdir()
        REPORT.mkdir()
    x, y = old_wind("zonal"), old_wind("merid")
    old_magnitude = np.array([np.hypot(x[i], y[i]).max() for i in range(768)])
    peak = int(old_magnitude.argmax())
    old_peak = float(old_magnitude[peak])
    spatial_x = np.asarray(x[peak]).copy() * (TAU / old_peak)
    spatial_y = np.asarray(y[peak]).copy() * (TAU / old_peak)
    # Verify that all original records factor into this same vector map.
    spatial_error = 0.
    for i in range(768):
        factor = old_magnitude[i] / TAU
        spatial_error = max(spatial_error, float(np.max(np.abs(x[i] - spatial_x * factor))),
                            float(np.max(np.abs(y[i] - spatial_y * factor))))
    assert spatial_error < 1e-14, spatial_error
    times = (np.arange(NREC) + 0.5) * SPACING
    envelope = gaussian(times)
    # At startup MITgcm brackets t=0 with the last and first periodic records.
    # Supply the negative-time Gaussian in the last slot, never reached again
    # during this 30-day experiment (31-day cycle).
    envelope[-1] = gaussian(-SPACING / 2)
    if args.review_only:
        write_review(envelope, old_magnitude, old_peak, spatial_x, spatial_y, spatial_error)
        return
    for component, spatial in (("zonal", spatial_x), ("merid", spatial_y)):
        path = NEW_INPUT / wind_name(component)
        with path.open("xb") as f:
            for value in envelope:
                np.asarray(spatial * value, dtype=">f8").tofile(f)
        assert path.stat().st_size == NREC * NX * NY * 8
        print(f"Wrote {path.name}", flush=True)
    obcs = NEW_INPUT / "obcs_31day"
    obcs.mkdir()
    for name in ("OBWt", "OBWs", "OBSt", "OBSs", "OBNt", "OBNs"):
        src = INPUT / "obcs_ts_from_initial_hydrography_16day" / (name + ".bin")
        size = src.stat().st_size // 768
        with src.open("rb") as f:
            first = f.read(size)
            for _ in range(767):
                assert f.read(size) == first, f"Nonconstant boundary record: {src}"
            assert not f.read(1)
        with (obcs / (name + ".bin")).open("xb") as f:
            for _ in range(NREC):
                f.write(first)
    diffs = []
    for base in ("Rplus", "Rminus"):
        dst = ROOT / (base + "_G5d")
        dst.mkdir()
        for name in FILES:
            src = ROOT / base / name
            if name not in ("data", "data.obcs"):
                shutil.copy2(src, dst / name)
                assert src.read_bytes() == (dst / name).read_bytes()
                continue
            old = src.read_text()
            new = old
            if name == "data":
                new = replace(new, "endTime", f"{END}.")
                new = replace(new, "externForcingCycle", f"{CYCLE}.")
                new = replace(new, "zonalWindFile", f"'../gaussian5d_input/{wind_name('zonal')}'")
                new = replace(new, "meridWindFile", f"'../gaussian5d_input/{wind_name('merid')}'")
            else:
                new = new.replace("../west_south_input/obcs_ts_from_initial_hydrography_16day/",
                                  "../gaussian5d_input/obcs_31day/")
            (dst / name).write_text(new)
            diffs.extend(difflib.unified_diff(old.splitlines(True), new.splitlines(True),
                                            fromfile=f"{base}/{name}", tofile=f"{dst.name}/{name}"))
        link_new(dst / "mitgcmuv", ROOT / "build" / "mitgcmuv")
    (REPORT / "configuration.diff").write_text("".join(diffs))
    write_review(envelope, old_magnitude, old_peak, spatial_x, spatial_y, spatial_error)


def write_review(envelope, old_magnitude, old_peak, spatial_x, spatial_y, spatial_error):
    t = np.arange(0, END, DT, dtype=float)
    applied = TAU * interpolate(envelope, t)
    # The original experiment ends before its 16-day cycle wraps. For a fair
    # longer comparison include its single pulse, not artificial repetitions.
    original_records = np.zeros(NREC)
    original_records[:len(old_magnitude)] = old_magnitude
    original_actual = interpolate(original_records, t)
    original_equal_peak = original_actual * TAU / old_peak
    analytic = TAU * gaussian(t)
    # Including every interpolation knot and both endpoints makes the
    # trapezoidal integral exact for the piecewise-linear applied forcing.
    knots = np.concatenate(([0.], (np.arange(END // SPACING) + 0.5) * SPACING, [float(END)]))
    def impulse(records):
        values = interpolate(records, knots)
        return float(np.dot(np.diff(knots), (values[:-1] + values[1:]) / 2))

    impulses = {
        "original_binary_actual": impulse(original_records),
        "original_equal_peak_0p1": impulse(original_records) * TAU / old_peak,
        "gaussian5d_applied": TAU * impulse(envelope),
        "gaussian5d_formula": TAU * SIGMA * math.sqrt(math.pi / 2) * (
            math.erf((END - TC) / (math.sqrt(2) * SIGMA))
            - math.erf(-TC / (math.sqrt(2) * SIGMA))),
    }
    iy, ix = np.unravel_index(np.hypot(spatial_x, spatial_y).argmax(), spatial_x.shape)
    direction_x, direction_y = spatial_x[iy, ix] / TAU, spatial_y[iy, ix] / TAU
    with (REPORT / "wind_stress_impulse.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["forcing", "start_s", "end_s", "magnitude_impulse_N_s_m2",
                         "magnitude_impulse_N_day_m2", "zonal_impulse_N_s_m2", "meridional_impulse_N_s_m2"])
        for name, value in impulses.items():
            writer.writerow([name, 0, END, value, value / 86400, value * direction_x, value * direction_y])
    # Unwindowed, non-demeaned finite-pulse transform. Common dt and duration;
    # no per-spectrum normalization or zero padding. Units: (N/m^2)^2 day^2.
    frequency = np.fft.rfftfreq(len(t), DT / 86400)
    power_old = np.abs(np.fft.rfft(original_equal_peak) * DT / 86400) ** 2
    power_new = np.abs(np.fft.rfft(applied) * DT / 86400) ** 2
    np.savetxt(REPORT / "forcing_spectrum.csv", np.column_stack((frequency, power_old, power_new)),
               delimiter=",", header="frequency_cpd,original_equal_peak_power,gaussian_applied_power", comments="")
    np.savetxt(REPORT / "forcing_timeseries.csv", np.column_stack((t, t / 86400, original_actual,
               original_equal_peak, analytic, applied)), delimiter=",", comments="",
               header="time_s,time_days,original_actual_N_m2,original_equal_peak_N_m2,gaussian_formula_N_m2,gaussian_applied_N_m2")
    key_t = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 10, 15, 30]) * 86400
    np.savetxt(REPORT / "gaussian_key_times.csv", np.column_stack((key_t / 86400,
               TAU * gaussian(key_t), TAU * interpolate(envelope, key_t))), delimiter=",", comments="",
               header="time_days,gaussian_formula_N_m2,gaussian_applied_N_m2")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, limit, unit in ((axes[0, 0], 8, 86400), (axes[0, 1], 8, 3600)):
        select = t <= limit * unit
        ax.plot(t[select] / unit, original_equal_peak[select], label="Original, scaled to common peak", color="C1")
        ax.plot(t[select] / unit, original_actual[select], label="Original binary (actual)", color="C1", ls=":")
        ax.plot(t[select] / unit, analytic[select], label="5-day Gaussian formula", color="C0")
        ax.plot(t[select] / unit, applied[select], label="5-day MITgcm interpolation", color="C2", ls="--")
        ax.set(xlabel="Time (days)" if unit == 86400 else "Time (hours)", ylabel="Wind stress (N/m$^2$)", xlim=(0, limit))
        ax.grid(alpha=.2)
    axes[0, 0].axvline(2.5, color=".6", lw=.7)
    axes[0, 0].axvline(5, color=".6", lw=.7, ls=":")
    axes[0, 0].legend(fontsize=8)
    ax = axes[1, 0]
    select = (frequency > 0) & (frequency <= 5)
    ax.semilogy(frequency[select], power_old[select], label="Original, common peak", color="C1")
    ax.semilogy(frequency[select], power_new[select], label="5-day applied Gaussian", color="C0")
    ax.axvline(1, color=".5", ls="--", label="1 cpd")
    ax.set(xlabel="Frequency (cycles/day)", ylabel="Squared Fourier amplitude ((N/m$^2$)$^2$ day$^2$)")
    ax.legend(fontsize=8)
    ax.grid(alpha=.2)
    magnitude = np.hypot(spatial_x, spatial_y)
    ax = axes[1, 1]
    im = ax.pcolormesh(magnitude, shading="auto", cmap="viridis", vmin=0, vmax=TAU)
    stride = 40
    yy, xx = np.mgrid[0:NY:stride, 0:NX:stride]
    ax.quiver(xx, yy, spatial_x[::stride, ::stride], spatial_y[::stride, ::stride])
    ax.set(xlabel="Grid column", ylabel="Grid row", title="Preserved vector map and buffer taper at envelope peak")
    fig.colorbar(im, ax=ax, label="N/m$^2$")
    fig.savefig(REPORT / "forcing_comparison.png", dpi=180)
    fig.savefig(REPORT / "forcing_comparison.pdf")
    plt.close(fig)
    inertial = int(np.argmin(np.abs(frequency - 1)))
    band = (frequency >= .8) & (frequency <= 1.2)
    metrics = dict(run_days=30, frequency_resolution_cpd=1/30, tc_seconds=TC, sigma_seconds=SIGMA,
                   original_binary_peak_N_m2=old_peak, new_formula_peak_N_m2=TAU,
                   new_applied_peak_N_m2=float(applied.max()),
                   power_ratio_at_1_cpd=float(power_new[inertial] / power_old[inertial]),
                   power_ratio_0p8_to_1p2_cpd=float(power_new[band].sum() / power_old[band].sum()),
                   vector_factorization_error_N_m2=spatial_error,
                   max_interpolation_error_first_10d_N_m2=float(np.max(np.abs(applied[:14400]-analytic[:14400]))))
    metrics["wind_stress_impulse_N_s_m2"] = impulses
    metrics["impulse_ratio_gaussian_to_original_equal_peak"] = impulses["gaussian5d_applied"] / impulses["original_equal_peak_0p1"]
    metrics["impulse_ratio_gaussian_to_original_actual"] = impulses["gaussian5d_applied"] / impulses["original_binary_actual"]
    (REPORT / "verification.json").write_text(json.dumps(metrics, indent=2) + "\n")
    readme = f"""# Gaussian wind sensitivity experiments

Rplus_G5d and Rminus_G5d copy Rplus and Rminus respectively, using the same
bathymetry, grid, hydrography, initial state, OBCS physics, mixing, executable,
60-s timestep, and 1200-s diagnostics. No production files were changed.

Wind envelope: exp[-(t-216000)^2/(2*72000^2)], nominal peak 0.1 N/m^2.
The original binary vector map (including buffer taper) is factored from its
largest nonzero record; only its scalar temporal multiplier changes.
The actual original sampled peak is {old_peak:.10f} N/m^2. The review includes
both this actual series and its rescaling to the common 0.1 N/m^2 peak.

Run duration: 30 days (2592000 s); frequency spacing: 0.033333 cpd.
Forcing record spacing stays 1800 s. MITgcm uses midpoint records and linear
interpolation, so the Gaussian's two samples bracketing day 2.5 give a tiny
plateau from day 2.489583 to 2.510417, peaking at {applied.max():.10f} N/m^2
(0.0078% below the analytic maximum). No cutoff is applied at day 2.5 or 5.
At day 5 the analytic stress is {TAU*gaussian(5*86400):.10f} N/m^2 and continues decaying.

The 31-day cycle prevents end-of-run wraparound. Its last record is a startup
guard evaluated at -900 s, used to bracket t=0; the run stops before this guard
is reached again. All other records sample the full Gaussian at midpoint times.
This yields a continuous interpolated startup and avoids a half-amplitude
artifact at t=0. Floating-point underflow occurs only for negligible late tails.

Boundary T/S files repeat the exact production records to 1488 records, matching
the shared forcing cycle. Only their paths and record counts change, not values.
data.pkg, data.diagnostics and eedata are byte-identical copies of the templates.
configuration.diff shows every changed namelist line.

Spectra use the same 30-day, 60-s sampled interval, no taper, no mean subtraction,
no zero padding, and squared Fourier amplitudes with time in days. The original
comparison pulse is interpolated from the actual binary record amplitudes and
scaled to the common peak, then continued at zero without repeating the original
16-day cycle. The Gaussian spectrum uses MITgcm interpolation.
New/original power at 1 cpd: {metrics['power_ratio_at_1_cpd']:.6g}.
Power ratio summed over 0.8--1.2 cpd: {metrics['power_ratio_0p8_to_1p2_cpd']:.6g}.
The nonzero starting Gaussian (1.1% of peak) contributes a small high-frequency
tail in a finite run; it is not the spectrum of an infinite Gaussian.
The longer pulse also supplies a much larger total wind impulse: this tests
duration/shape at fixed peak, not fixed integrated stress.
The scientific comparison targets the dominant frequency of the bay-generated
response rather than its absolute amplitude. No impulse normalization is applied
to either simulation; the configured Gaussian and run inputs remain unchanged.

Time-integrated stress from t=0 to day 30, at an untapered active-domain point:

| Forcing | Impulse (N s/m^2) | Impulse (N day/m^2) |
| --- | ---: | ---: |
| Original binary, actual | {impulses['original_binary_actual']:.6f} | {impulses['original_binary_actual']/86400:.9f} |
| Original, common 0.1 peak | {impulses['original_equal_peak_0p1']:.6f} | {impulses['original_equal_peak_0p1']/86400:.9f} |
| 5-day Gaussian, applied | {impulses['gaussian5d_applied']:.6f} | {impulses['gaussian5d_applied']/86400:.9f} |
| 5-day Gaussian, formula | {impulses['gaussian5d_formula']:.6f} | {impulses['gaussian5d_formula']/86400:.9f} |

The applied Gaussian supplies {metrics['impulse_ratio_gaussian_to_original_equal_peak']:.6f}
times the impulse of the common-peak original pulse, or
{metrics['impulse_ratio_gaussian_to_original_actual']:.6f} times that of the actual
production binary pulse. Applied integrals use exact trapezoidal integration
of all MITgcm interpolation knots and both endpoints. The formula integral uses
the Gaussian error function. The original pulse is continued at zero after its
decay, without periodic repetition. These are stress impulses per unit area,
not domain-integrated forces. In the buffers, multiply by the unchanged spatial
taper. wind_stress_impulse.csv also records signed zonal/meridional components
at the untapered point; verification.json stores the scalar integrals and ratios.

Plots: forcing_comparison.png/pdf. Data: forcing_timeseries.csv,
gaussian_key_times.csv, forcing_spectrum.csv, wind_stress_impulse.csv.
Checks: verification.json.

MPI decomposition requires exactly 16 ranks (4 x 4), within the 20-core limit.
The launcher runs only the two new cases and refuses directories with outputs.
Output is written to LaCie. Existing diagnostics imply about 2.58 TB per case
for 30 days, excluding dumps/checkpoints; allow approximately 5.3 TB for the pair.
LaCie has approximately 5.6 TB available at setup. A 40-day pair would need
6.88 TB for diagnostics alone, so 30 days was selected to fit the drive.
Runtime is approximately twice a 15-day case at the same throughput.

Launch after reviewing the figures:

```bash
bash CTW_SVB/scripts/run_west64_south72_g5d_30d_sequential.sh
```

No run has been launched by this setup.
"""
    (REPORT / "README.md").write_text(readme)
    print(json.dumps(metrics, indent=2), flush=True)


if __name__ == "__main__":
    main()
