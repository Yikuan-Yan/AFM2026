# AFM2026

Reproducible analysis code, derived results, quality-control figures and experimental planning for silica colloidal-probe AFM measurements in water–glycerol mixtures at 25.6 °C.

## Repository contents

- `analysis/calibrate_cantilevers.py`: raw thermal-noise and hard-contact calibration workflow.
- `analysis/fit_glycerol_surface_forces.py`: raw reconstruction and silica sphere–plane surface-force fitting.
- `analysis/analyze_velocity_systematics.py`: contact, snap-in, pull-off, raster and acquisition-history diagnostics.
- `analysis/analyze_velocity_joint_fit.py`: same-pixel velocity comparisons, hydrodynamic shape tests and joint PB fits.
- `analysis/*_results/`: derived CSV data, figures, reports, provenance and SHA-256 manifests.
- `analysis/FOLLOWUP_8X8_PALINDROME_FORCE_MAPPING_PROTOCOL.md`: detailed follow-up experiment using balanced 8×8 palindromic force mapping.

The current scientific interpretation and its limitations are documented in:

- [`analysis/velocity_joint_fit_results/REPORT.md`](analysis/velocity_joint_fit_results/REPORT.md)
- [`analysis/velocity_systematics_results/REPORT.md`](analysis/velocity_systematics_results/REPORT.md)
- [`analysis/surface_force_results/REPORT.md`](analysis/surface_force_results/REPORT.md)
- [`analysis/results/REPORT.md`](analysis/results/REPORT.md)

## Raw-data policy

Original JPK force maps, individual force curves, thermal-noise files and acquisition photographs are intentionally not stored in this Git repository. The `.gitignore` excludes the local `20-08-26/` and `calibration/` trees as well as the raw JPK/TND file extensions globally.

The derived provenance and input manifests retain source identities and hashes so that results can be checked against an authorized raw-data archive without placing the raw measurements in Git.

## Reproduction

With the excluded raw folders restored at the paths expected by the scripts, the analyses are run from the repository root:

```bash
python analysis/calibrate_cantilevers.py
python analysis/fit_glycerol_surface_forces.py
python analysis/analyze_velocity_systematics.py
python analysis/analyze_velocity_joint_fit.py
```

Each result directory contains an `artifact_manifest.sha256` file for checking the committed derived artifacts.
