#!/usr/bin/env bash
set -euo pipefail

run_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$run_dir"

for required_file in long_term_tracking.init lat.bmad fccee_z.bmad laser_ramper_expression.bmad; do
  if [[ ! -s "$required_file" ]]; then
    echo "ERROR: missing required file: $required_file" >&2
    exit 2
  fi
done

# The DM module looks for these three table files in the working directory
# by default. Symlink them in from the canonical copies committed under
# bmad/laser/ on this branch, rather than duplicating ~21 MB of table data
# per run directory.
dm_laser_dir="$run_dir/../../bmad/laser"
for table in dm_energy_table.dat dm_spin_table.dat dm_cdf_table.dat; do
  if [[ ! -s "$table" ]]; then
    if [[ -s "$dm_laser_dir/$table" ]]; then
      ln -s "$dm_laser_dir/$table" "$table"
    else
      echo "ERROR: $table not found here or at $dm_laser_dir/$table" >&2
      echo "       (expected to find it committed under bmad/laser/ on the density-matrix branch)" >&2
      exit 2
    fi
  fi
done

if [[ -n "${BMAD_ENV_SCRIPT:-}" ]]; then
  # shellcheck disable=SC1090
  source "$BMAD_ENV_SCRIPT"
fi

ltt_exe=${LTT_EXE:-long_term_tracking}
if ! command -v "$ltt_exe" >/dev/null 2>&1 && [[ ! -x "$ltt_exe" ]]; then
  echo "ERROR: long_term_tracking not found. Set LTT_EXE to the custom executable." >&2
  exit 3
fi

echo "Run directory: $run_dir"
echo "Executable: $(command -v "$ltt_exe" 2>/dev/null || printf '%s' "$ltt_exe")"
echo "Input SHA256: $(shasum -a 256 long_term_tracking.init | awk '{print $1}')"
"$ltt_exe" long_term_tracking.init
