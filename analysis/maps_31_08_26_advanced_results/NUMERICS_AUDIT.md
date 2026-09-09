# Numerical audit: 31-08-26 advanced analysis

## Computational claim

- Claimed output: model-light same-pixel shape/spatial evidence plus apparent sphere-plane PB parameters.
- Domain: four same-location 8x8 maps, 256 curves, 20-200 nm in 5 nm bins, 25.6 C pure-water environment.
- Units: distance nm, force pN, potential mV, time min, k N/m, InvOLS nm/V.
- Numerical contract: formula/prediction closure <=1e-7 pN; qualitative inference uses spatial and model sensitivity.
- Reproducibility: deterministic except seeded spatial permutation/bootstrap; seed=20260901.

## Direct probes

- Nonlinear-PB library identity against the existing sphere-plane evaluator: max `9.095e-13 pN` — PASS.
- Nonlinear PB to linear-HHF limit at `|psi|=0.1 mV`, `lambda=20 nm`: max relative difference `1.656e-03` — PASS.
- Fixed lubrication limiting laws: relative range of `F*D` `1.656e-16`; max relative error in `F(2v)=2F(v)` `0.000e+00` — PASS.
- Saved prediction/residual algebraic closure: max `0.000e+00 pN` — PASS.
- Endpoint normalization closure, valid curves: max `0.000e+00` — PASS.
- Independent-fit RSS is no larger than every parameter-sharing model: `True` — PASS.
- All 256 primary fitted parameter and diagnostic arrays finite: `True` — PASS.

## Physical-scale checks

At median measured gap speed `1.995909 µm/s`, the fixed ideal no-slip sphere-plane lubrication term is 20 nm: 34.25 pN, 50 nm: 13.70 pN, 100 nm: 6.85 pN, 200 nm: 3.42 pN.
The fitted dimensionless potential `e|psi|/kBT` has median `0.913` [IQR `0.680`, `1.146`]. The primary evaluator is nonlinear PB, so this is a regime descriptor rather than a linearization gate.

## Claim limits

- Adjacent 5 nm bins from one curve are correlated; profile regions and AIC/BIC are diagnostic, not formal independent-observation confidence statements.
- The nonlinear model assumes equal constant potential, symmetric 1:1 PB electrolyte, ideal Derjaguin geometry, fixed Hamaker constant, and ideal no-slip hydrodynamics.
- Spring-constant common scaling, contact-zero shifts, baseline choice, spatial clustering, and model form are audited separately; no single branch is promoted to a bulk Debye-length or zeta-potential measurement.
- InvOLS uncertainty changes both force and separation and is not represented by the common-force-scale branch alone.
