# Gauss_Sim

Calculator and simulator for the curvature of a surface. The input is only lengths measured on the surface. The program finds the Gaussian curvature K and the radius R. No angle and no radius is an input.

## Build and run

```
cargo run --release
cargo test
```

## Method: cut the triangle, cancel the angle

1. Measure the three sides a, b, c of a triangle ABC on the surface. Three lengths alone do not fix K. A triangle with these sides exists on every surface of constant K (up to a size limit for K > 0).
2. Put a mark D on one side. Measure m along that side from its start Q. Measure the cut d from the opposite vertex P to D.
3. The cut makes a sub-triangle PQD. It shares the angle at Q with the full triangle PQQ′. Write the law of cosines (1) for both triangles and eliminate that angle. The result (2) holds only lengths and K.
3x. Cross cut: join a mark X on PQ to a mark Y on PQ′. The small triangle PXY shares the angle at P with PQQ′. Eliminating it gives (2x). Marks at the midpoints of both sides give d₀ = ℓ/2 on a plane.
4. Solve (2) or (2x) for K. Each cut gives one value Kⱼ. More cuts (other positions, other apexes) give more values, combined into K*. If they agree (χ²/dof ≈ 1), one curvature fits the whole triangle.
5. On a plane, (2) becomes Stewart's theorem (3). The ratio ρ = d/d₀ is the direct test: ρ = 1 is flat, ρ > 1 is sphere-like, ρ < 1 is saddle-like.
6. From K*: R = 1/√|K*|, and ω = 1/R is the turn in radians per unit distance. The angles α, β, γ and the angle excess ε are outputs.

Gauss (1828, Art. 12) showed that K depends only on lengths measured on the surface. This program uses only such lengths.

## Notation

| Symbol | Meaning |
|---|---|
| A, B, C | vertices of the triangle (marks on the surface) |
| a, b, c | measured sides: a = \|BC\|, b = \|CA\|, c = \|AB\| |
| Pⱼ | apex of cut j (A, B or C) |
| (P, Q, Q′) | cyclic order: (A, B, C), (B, C, A), (C, A, B). The base of cut j is QQ′, the side opposite P |
| foot cut | from the apex P to a mark Dⱼ on the base QQ′ |
| cross cut | from a mark X on PQ to a mark Y on PQ′ (both sides that meet at P) |
| Dⱼ | foot of cut j, a mark on the base |
| ℓⱼ = \|QQ′\| | base length (a, b or c) |
| uⱼ = \|PQ\|, vⱼ = \|PQ′\| | the two other sides, from the apex |
| mⱼ = \|QDⱼ\| | measured along the base, from Q |
| nⱼ = ℓⱼ − mⱼ | rest of the base, \|DⱼQ′\| |
| dⱼ = \|PⱼDⱼ\| | measured cut length |
| xⱼ = \|PX\|, yⱼ = \|PY\| | cross cut: marks measured from P along PQ and PQ′ |
| dⱼ = \|XY\| | cross cut length |
| tⱼ | position: m/ℓ (foot) or x/u (cross) |
| σ | standard deviation of one length measurement (same for every length) |
| K | Gaussian curvature, unit⁻². K > 0 sphere-like, K = 0 flat, K < 0 saddle-like |
| Kⱼ | curvature from cut j alone |
| K* | combined curvature: inverse-variance mean of the Kⱼ |
| σ(Kⱼ), σ(K*) | standard deviations from σ; σ(K*) uses every length once, so shared sides count once |
| K_sim, R_sim, a_sim, b_sim, c_sim | simulator surface and true sides (before noise) |
| R = 1/√\|K*\| | radius of curvature (sphere radius, or saddle pseudo-radius) |
| Rⱼ = 1/√\|Kⱼ\| | radius from cut j alone |
| ω = √K* = 1/R | turn per unit distance, rad/unit (K* > 0) |
| subscript ₖ | evaluated at curvature K: snₖ, csₖ, Gₖ |
| snₖ(x), csₖ(x) | sin(√K x)/√K, cos(√K x) for K > 0; x, 1 for K = 0; sinh(√−K x)/√−K, cosh(√−K x) for K < 0 |
| d₀ | cut length that a plane (K = 0) predicts from the other lengths |
| Δ = d² − d₀² | defect, unit². Zero on a plane |
| ρ = d / d₀ | ratio test. ρ = 1 flat, ρ > 1 for K > 0, ρ < 1 for K < 0 |
| d(K) | cut length that (2) predicts on curvature K; d₀ = d(0) |
| Gₖ | cut function; its root is Kⱼ |
| α, β, γ | angles at A, B, C. Outputs only, computed from lengths and K* |
| ε = α + β + γ − π | angle excess, ε = K𝒜 |
| 𝒜 | triangle area, unit² |

## Equations

```
(1) law of cosines   csₖ(a) = csₖ(b)·csₖ(c) + K·snₖ(b)·snₖ(c)·cos α
(2) cut relation     csₖ(v)·snₖ(m) + csₖ(u)·snₖ(n) = csₖ(d)·snₖ(ℓ)
(2x) cross relation  [csₖ(ℓ) − csₖ(u)csₖ(v)] / [K snₖ(u)snₖ(v)] = [csₖ(d) − csₖ(x)csₖ(y)] / [K snₖ(x)snₖ(y)]
(3) flat limit       v²m + u²n = ℓ(d² + mn)   ⇒   d₀² = (v²m + u²n)/ℓ − mn
(3x) flat cross cut  d₀² = x² + y² − xy(u² + v² − ℓ²)/(uv)
(4) defect           Δ = d² − d₀²
(5) cut function     Gₖ = [csₖ(v)snₖ(m) + csₖ(u)snₖ(n) − csₖ(d)snₖ(ℓ)] / K,   G₀ = ℓΔ/2
(5x) cross function  Gₖ = left side of (2x) − right side of (2x)
(6) solve            Gₖ(Kⱼ) = 0,   sign Kⱼ = sign Δⱼ = sign(ρⱼ − 1)
(7) Gauss–Bonnet     ε = α + β + γ − π = K𝒜
```

Both sides of (2x) equal cos of the angle at P (law of cosines (1) for PQQ′ and for PXY). The program writes each side in a form exact at K = 0.

Derivation of (2): the law of cosines (1) at Q, for triangle PQQ′ and for triangle PQD, gives two expressions for K·snₖ(u)·cos Q. Set them equal. The term csₖ(m)snₖ(ℓ) − csₖ(ℓ)snₖ(m) is snₖ(ℓ − m) = snₖ(n).

The program evaluates (5) with the Stumpff functions c₂(z) = (1 − cos√z)/z and c₃(z) = (√z − sin√z)/z^{3/2}, with z = Kx²:

```
csₖ(x) = 1 − z·c₂(z),   snₖ(x) = x·(1 − z·c₃(z))
Gₖ = ℓ³c₃(Kℓ²) − m³c₃(Km²) − n³c₃(Kn²) − v²c₂(Kv²)snₖ(m) − u²c₂(Ku²)snₖ(n) + d²c₂(Kd²)snₖ(ℓ)
```

This form has no cancellation at K = 0, so it holds for survey triangles where K·ℓ² is near 10⁻⁴ or less.

The cut length on a surface of curvature K increases with K. This is the triangle comparison of Toponogov (1959) and Alexandrov (see Burago, Burago & Ivanov, 2001; Bridson & Haefliger, 1999). So Gₖ > 0 below the root and Gₖ < 0 above it, the root is unique, and bisection finds it. For K > 0 the search stops below (2π/perimeter)².

## Uncertainty and verdict

- σ(Kⱼ) comes from σ on each of the cut's lengths ℓ, u, v, m (or x, y), d (central differences).
- K* is the inverse-variance mean of the Kⱼ.
- σ(K*) treats K* as one function of every measured length (a, b, c and each mⱼ, dⱼ). The cuts share the sides, so this counts each side once. On flat noisy data the 2σ verdict reads "curved" about 5% of the time, as expected.
- χ²/dof tests if the Kⱼ agree. Shared sides make the Kⱼ errors correlate, so χ²/dof is a little below 1 on good data.
- Verdict: |K*| ≤ 2σ(K*) is flat within measurement. Then the program gives the smallest sphere radius and saddle pseudo-radius that the data allow at 2σ(K*).
- With exact data (σ = 0), σ(K*) is the scatter of the Kⱼ, and |K*|·L² below 10⁻¹² counts as zero (floating-point floor).
- The solver searches K·L² ≥ −400 (R > L/20, L the largest length). Below that the cut function loses its sign in double precision.

## Tabs

- **Results**: verdict, K* ± σ(K*), R, ω, χ²/dof, derived angles, excess, area, per-cut table.
- **Ratio test ρ**: ρ − 1 for each cut against foot position t, with ±2σ bars and the model curve for K*.
- **Cut functions Gₖ**: Gₖ(K) for each cut. Each curve crosses zero once.
- **Intrinsic map**: azimuthal-equidistant map about A. Distances from A and angles at A are true. The gap between side BC and the dashed chord shows curvature. A display factor can enlarge K* for the drawing only.
- **Sphere view**: the triangle on a sphere of radius R (K* > 0). Drag to rotate, scroll to zoom. A complete surface with K < 0 has no smooth isometric embedding in 3-D space (Hilbert, 1901).
- **Notation & method**: the tables above.
- **Report**: plain-text execution trace: inputs, assumptions, each cut's flat prediction with the arithmetic, the full bisection log (bracket, K tried, Gₖ), every snₖ/csₖ term of the equation at the solved K, the combination sums, and the derivation of R, ω, the angles and the area. It ends with a standalone Python script (mpmath, 50 digits, plain sin/cos/sinh/cosh, none of this program's code) that recomputes every Kⱼ, K*, R and the angles. "Copy report" puts it on the clipboard; "Save report as .txt" writes `~/Documents/gauss_sim_report_<time>.txt`.

Measurements are saved to `~/.local/share/gauss_sim/measurements.txt` on every edit and loaded at start.

## Presets

- 100 km survey triangle on R = 6371 km, σ = 1 cm.
- Gauss (1828, Art. 28) triangle: A Hohehagen, B Brocken, C Inselberg. Sides about 107, 85 and 69 km (translators' note in Gauss, 1902), with Brocken–Inselberg the longest side: a ≈ 107, b ≈ 85, c ≈ 69 km. Gauss gives an angle excess of 14.85348″. With these rounded sides on R = 6371 km the excess is 14.887″.
- Octant of a sphere, R = 1, exact data: every cut gives K = 1.
- Flat 1 km field, σ = 1 mm.
- Saddle with pseudo-radius 10.
- 1 km triangle on R = 6371 km, σ = 1 mm: flat within measurement.

## Sources

- Bridson, M. R., & Haefliger, A. (1999). *Metric spaces of non-positive curvature* (Grundlehren der mathematischen Wissenschaften, Vol. 319). Springer. https://doi.org/10.1007/978-3-662-12494-9
- Burago, D., Burago, Y., & Ivanov, S. (2001). *A course in metric geometry* (Graduate Studies in Mathematics, Vol. 33). American Mathematical Society. https://doi.org/10.1090/gsm/033
- Danby, J. M. A. (1988). *Fundamentals of celestial mechanics* (2nd ed.). Willmann-Bell.
- Gauss, C. F. (1828). Disquisitiones generales circa superficies curvas. *Commentationes Societatis Regiae Scientiarum Gottingensis Recentiores*, 6, 99–146.
- Gauss, K. F. (1902). *General investigations of curved surfaces of 1827 and 1825* (J. C. Morehead & A. M. Hiltebeitel, Trans.). Princeton University Library.
- Hilbert, D. (1901). Ueber Flächen von constanter Gaussscher Krümmung. *Transactions of the American Mathematical Society*, 2(1), 87–99. https://doi.org/10.1090/S0002-9947-1901-1500557-5
- Stewart, M. (1746). *Some general theorems of considerable use in the higher parts of mathematics*. Sands, Murray & Cochran.
- Stumpff, K. (1959). *Himmelsmechanik* (Vol. 1). VEB Deutscher Verlag der Wissenschaften.
- Toponogov, V. A. (1959). Riemannian spaces with curvature bounded below. *Uspekhi Matematicheskikh Nauk*, 14(1), 87–130.
