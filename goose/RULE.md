# Pre-declared analysis rule (written 2026-10-23 10:55 MDT, before the run)

Quantity: K_intrinsic from  eps_obs = sum(observed angles) - 180 deg  =  K * A.

Area A for each triangle comes from PRINTED SIDES (ground measurements), not coordinates:
  3 sides -> Heron;  1-2 sides -> A = s^2 sin(a1) sin(a2) / (2 sin(a3)), a3 = angle opposite s.
  Units: feet = 0.3048 m, toise = 1.949036 m, logs as printed (US 39th: log metres; Helmert: log metres;
  Rhodesia: log English feet; Anglo-French: log toises). Triangles whose side-area disagrees with
  coordinate-area by more than 2% are flagged as transcription/unit faults and fall to tier C.

Provenance tiers (per triangle; weight multiplier):
  A (1.00)  three observed angles + all three printed sides  -> rebuildable from scratch, overdetermined
  B (0.50)  three observed angles + one or two printed sides
  C (0.25)  three observed angles, area from coordinates only
  excluded  angles computed/adjusted (Brousseaud "spherical"), missing, or inconsistent with printed sum

Noise: sigma_s per survey = 1.4826 * MAD of residuals from that survey's own one-parameter fit eps = k_s A
       (no ellipsoid used).
Signal-to-noise: SNR_i = K0 * A_i / sigma_s, K0 = 1/(6371 km)^2 (a fixed scale, not a fit).
Cut: keep SNR_i >= 2.
Fit: weighted least squares through origin, eps = K A, weights w_i = tier_i / sigma_s^2.
Error: sandwich variance for non-optimal weights, inflated by sqrt(chi2/dof) if > 1.
Headline = the run above. Sensitivity runs (reported, not chosen between): SNR cut 0/1/2/3;
tier weights equal / A-only; area from coordinates instead of sides.
Comparison model: K_ell = 1/(MN) on WGS84 at each triangle's mean latitude.
