# Written by Claude, approved by KRM.
"""Delete the CTW_SVB output approved on 2026-10-29.

Dry run by default (lists counts and sizes).  Pass --execute to delete.

  LaCie  : timestamped *.data/*.meta of R_G5dplus and R50Lplus (logs, grid,
           pickups kept); folders smoke_Rplus_50L_hfac, Splus, Sminus.
  local  : smoke-test output (*.data/*.meta, mnc_*, per-rank STDOUT/STDERR);
           five superseded build folders; *.bin of the superseded-domain inputs.

Rplus/Rminus on LaCie, AmeliaProcessedOutput and anything outside the paths
below are never touched.
"""
import os
import re
import sys

L = "/Volumes/LaCie/CTW_SVB_reruns_telescoped_west64_south72"
C = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAMPED = re.compile(r"^[A-Za-z0-9_]+\.\d{10}\.(data|meta)$")
SMOKE = re.compile(r"(\.(data|meta|nc)$|^STDERR\.|^STDOUT\.(000[1-9]|00[1-9]\d)$)")

files, dirs = [], []
for run in ("R_G5dplus", "R50Lplus"):
    for e in os.scandir(os.path.join(L, run)):
        if e.is_file(follow_symlinks=False) and STAMPED.match(e.name):
            files.append(e.path)
dirs += [os.path.join(L, d) for d in ("smoke_Rplus_50L_hfac", "Splus", "Sminus")]

for d in ("archive/2_telescoped_buffer/ctw_reruns_telescoped/smoke_Rplus_binary",
          "archive/2_telescoped_buffer/ctw_reruns_telescoped/smoke_Rplus_noexch2",
          "archive/3_telescoped_south/ctw_reruns_telescoped_south/smoke_Rplus",
          "ctw_reruns_telescoped_west64_south72/smoke_Rplus"):
    for e in os.scandir(os.path.join(C, d)):
        if e.is_symlink():
            continue
        if e.is_file() and SMOKE.search(e.name):
            files.append(e.path)
        elif e.is_dir() and e.name.startswith("mnc_"):
            dirs.append(e.path)
dirs += [os.path.join(C, d) for d in (
    "archive/1_original_512x612/ctw_reruns/build",
    "archive/2_telescoped_buffer/ctw_reruns_telescoped/build_binary",
    "archive/2_telescoped_buffer/ctw_reruns_telescoped/build_east",
    "archive/2_telescoped_buffer/ctw_reruns_telescoped/build_noexch2",
    "archive/3_telescoped_south/ctw_reruns_telescoped_south/build_binary")]
for d in ("input/telescoped_buffer", "input/telescoped_south_extended"):
    for r, _, fs in os.walk(os.path.join(C, d)):
        files += [os.path.join(r, f) for f in fs if f.endswith(".bin")]

for p in files + dirs:
    assert p.startswith((L + "/", C + "/")), p
    assert not p.startswith((L + "/Rplus/", L + "/Rminus/")), p


def size(p):
    if os.path.islink(p):
        return 0
    if os.path.isfile(p):
        return os.stat(p).st_size
    return sum(os.lstat(os.path.join(r, f)).st_size for r, _, fs in os.walk(p) for f in fs)


lacie = sum(size(p) for p in files + dirs if p.startswith(L))
local = sum(size(p) for p in files + dirs if p.startswith(C))
print(f"{len(files)} files, {len(dirs)} folders: LaCie {lacie/1e12:.2f} TB, local {local/1e9:.1f} GB")
if "--execute" not in sys.argv:
    print("dry run; pass --execute to delete")
    sys.exit()
for p in files:
    os.remove(p)
    ad = os.path.join(os.path.dirname(p), "._" + os.path.basename(p))
    if os.path.lexists(ad):
        os.remove(ad)
def rmtree(top):
    # On LaCie (exFAT) macOS removes the "._*" companion of a file together
    # with the file, so entries can vanish mid-walk; skip those.
    for r, ds, fs in os.walk(top, topdown=False):
        for n in fs + [d for d in ds if os.path.islink(os.path.join(r, d))]:
            try:
                os.remove(os.path.join(r, n))
            except FileNotFoundError:
                pass
        for n in ds:
            p = os.path.join(r, n)
            if not os.path.islink(p) and os.path.isdir(p):
                os.rmdir(p)
    os.rmdir(top)


for d in dirs:
    if not os.path.lexists(d):
        continue
    rmtree(d)
    ad = os.path.join(os.path.dirname(d), "._" + os.path.basename(d))
    if os.path.lexists(ad):
        os.remove(ad)
print("deleted")
