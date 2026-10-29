#!/usr/bin/env python3
"""Create small MITgcm smoke-test run directories for CTW reruns."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
RERUNS = REPO / "CTW_SVB" / "ctw_reruns"
EXTERNAL_ROOT = Path("/Volumes/LaCie/CTW_SVB_reruns")


def replace_line(text: str, key: str, replacement: str) -> str:
    lines = []
    replaced = False
    for line in text.splitlines():
        if line.strip().startswith(key):
            lines.append(replacement)
            replaced = True
        else:
            lines.append(line)
    if not replaced:
        raise ValueError(f"Could not find line starting with {key!r}")
    return "\n".join(lines) + "\n"


def disable_diagnostics(text: str) -> str:
    for key in ("frequency(1)", "frequency(2)", "frequency(3)"):
        text = replace_line(text, key, f"  {key} = 0.,")
    return text


def create_case(name: str, wind: str, obcs_balance: str, delta_t: float, hours: float) -> Path:
    source = RERUNS / "Rplus"
    target = EXTERNAL_ROOT / name
    if target.exists():
        raise FileExistsError(f"{target} already exists; choose a new smoke case name")
    target.mkdir(parents=True)

    for item in ("data", "data.obcs", "data.diagnostics", "data.pkg", "data.mnc", "data.exch2", "eedata"):
        shutil.copy2(source / item, target / item)

    for link_name, link_target in (
        ("mitgcmuv", RERUNS / "build" / "mitgcmuv"),
        ("bathymetry", RERUNS / "bathymetry"),
        ("common_input", RERUNS / "common_input"),
    ):
        (target / link_name).symlink_to(link_target)

    data = (target / "data").read_text()
    data = replace_line(data, "endTime", f" endTime={hours * 3600:.0f}.,")
    data = replace_line(data, "deltaT", f" deltaT={delta_t:g}.,")
    data = replace_line(data, "chkptFreq", " chkptFreq=0.,")
    data = replace_line(data, "dumpFreq", " dumpFreq=0.,")
    data = data.replace(" monitorFreq=1200.,", " monitorFreq=1200.,\n dumpInitAndLast=.FALSE.,\n writePickupAtEnd=.FALSE.,")
    if wind == "off":
        data = replace_line(data, "meridWindFile", " meridWindFile=' ',")
        data = replace_line(data, "zonalWindFile", " zonalWindFile=' ',")
    (target / "data").write_text(data)

    obcs = (target / "data.obcs").read_text()
    if obcs_balance == "off":
        obcs = replace_line(obcs, "useOBCSbalance", " useOBCSbalance=.FALSE.,")
    (target / "data.obcs").write_text(obcs)

    diagnostics = disable_diagnostics((target / "data.diagnostics").read_text())
    (target / "data.diagnostics").write_text(diagnostics)
    return target


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("name")
    parser.add_argument("--wind", choices=("on", "off"), default="on")
    parser.add_argument("--obcs-balance", choices=("on", "off"), default="on")
    parser.add_argument("--delta-t", type=float, default=60.0)
    parser.add_argument("--hours", type=float, default=4.0)
    args = parser.parse_args()
    target = create_case(args.name, args.wind, args.obcs_balance, args.delta_t, args.hours)
    print(target)


if __name__ == "__main__":
    main()
