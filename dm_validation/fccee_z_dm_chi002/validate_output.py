#!/usr/bin/env python3
"""Validate the averages output from the DM chi_e=0.002 validation run.

Compares turn-0 vs turn-100000 live particle counts (the central question
this run exists to answer -- does DM stochastic recoil cause loss the LCFA
mean-kick run didn't show) and prints the final <Sy> for comparison against
reference/xip05_2_out_cons_LCFA_reference.ave's 0.006893103 (0.6893%).
"""

from pathlib import Path
import re
import sys


path = Path(sys.argv[1] if len(sys.argv) > 1 else "out_dm_chi002.ave")
if not path.is_file():
    raise SystemExit(f"ERROR: output not found: {path}")

header = {}
rows = []
for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
    if line.startswith("#"):
        match = re.match(r"#\s*([^=]+?)\s*=\s*(.*?)\s*$", line)
        if match:
            header[match.group(1).strip()] = match.group(2).strip()
        continue
    fields = line.split()
    if len(fields) >= 7:
        rows.append((int(fields[0]), int(fields[1]), float(fields[5])))

errors = []
if header.get("ltt%n_turns") != "100000":
    errors.append(f"ltt%n_turns is {header.get('ltt%n_turns')!r}, expected '100000'")
if header.get("ltt%rfcavity_on") != "T":
    errors.append("RF was not reported on")
if header.get("bmad_com%spin_tracking_on") != "T":
    errors.append("spin tracking was not reported on")
if not rows:
    errors.append("no averages rows found")
else:
    first, last = rows[0], rows[-1]
    if first[0] != 0:
        errors.append(f"first averages turn is {first[0]}, expected 0")
    if last[0] != 100000:
        errors.append(f"last averages turn is {last[0]}, expected 100000")
    print(f"turn-0 live particles: {first[1]} (requested: 2214)")
    print(f"turn-100000 live particles: {last[1]}")
    print(f"tracking losses after turn 0: {first[1] - last[1]}")
    print(f"final <Sy>: {last[2]:.9f} ({100*last[2]:.6f}%)")
    print("LCFA reference (different wavelength-equivalent point, n_particle=920): "
          "<Sy>=0.006893103 (0.6893%), 0 losses (920/920 retained)")
    if first[1] - last[1] > 0:
        print("*** DM run shows particle loss where the comparable LCFA mean-kick "
              "run showed none -- this is the effect this validation exists to detect. ***")

if errors:
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    raise SystemExit(1)
print("Basic output validation passed.")
