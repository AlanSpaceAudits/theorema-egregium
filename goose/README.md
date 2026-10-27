# Theorema Egregium: curvature measured from surveyors' triangles

Nine historical triangulation surveys (1739–1898, four continents) re-read from inside the surface only, as Gauss's Theorema Egregium allows. Each triangle's observed angle excess, divided by its area, gives the Gaussian curvature of the ground it was drawn on. Areas come from the printed side lengths, which are themselves carried from measured base lines. No latitudes, azimuths or stars go into the curvature result. The ellipsoid's curvature 1/(MN) is used only for comparison.

## 📄 [Read the report](https://rooster-ninja.github.io/TheoremaEgregium_CurvatureMeasuredFromSurveyorsTriangles/)

Hosted on GitHub Pages from `index.html`. It links directly to the spreadsheet in [`supplemental/`](https://github.com/rooster-ninja/TheoremaEgregium_CurvatureMeasuredFromSurveyorsTriangles/tree/main/supplemental) and the page extracts in [`sources/`](https://github.com/rooster-ninja/TheoremaEgregium_CurvatureMeasuredFromSurveyorsTriangles/tree/main/sources).

## Results

| Result | Value | Uncertainty |
| --- | --- | --- |
| Triangles in the result (signal ≥ 2× noise, rule fixed before the run) | 391 of 935 usable | — |
| Gaussian curvature K from observed angles | **2.450 × 10⁻¹⁴ m⁻²** | ± 0.012 (± 0.019 with survey-to-survey scatter) |
| Ellipsoid curvature 1/(MN), same triangles | 2.462 × 10⁻¹⁴ m⁻² | — |
| Ratio K / K ellipsoid | **0.995** | ± 0.005 (−1.0σ); ± 0.008 with survey scatter (−0.6σ) |
| Radius 1/√K | **6 389 km** | 1σ 6 373–6 404 km |
| Flat surface (K = 0) | rejected | > 100σ |
| US 39th parallel alone (84% of the weight) | ratio 1.000 | ± 0.005 |
| Whole-chain sums, core three surveys (373 triangles, 478 766 km²) | 2 374″ observed vs 2 432″ expected, ratio 0.977 | ± 0.011 |
| Constant per-triangle offset (fitted after the fact) | −0.15″ per triangle; ratio becomes 1.001 | ± 0.05″; ± 0.005 |

Within its error bars, the curvature measured only from angles and lengths on the ground agrees with the curvature of the ellipsoid the ground sits on. That is the empirical content of the Theorema Egregium. Latitude dependence of K (0.7% between 13° and 52°) is below the resolution. Two surveys read low on their own: the Eastern Oblique Arc (0.970 ± 0.015) and Helmert's 52° arc (0.958 ± 0.025). Appendix A and B suggest a small uniform offset in the printed angle sums (≈ −0.15″ per triangle) rather than a curvature difference. All results are fallible, inductive measurements, good to the sigma stated.

## Methods: claimed accuracy against what the data show

| Item | Accuracy claimed in the source | What the triangles' own closures show |
| --- | --- | --- |
| Epping base, 1857, 8 715.942 m (Bache–Würdemann bars) | p.e. ± 0.0158 m (1 in 552 000) | not testable from angles; enters only through area |
| Base apparatus (Hassler / Bache–Würdemann) | 1 : 200 000 / 1 : 310 000 | area error 10³–10⁵ × below the curvature signal |
| Single direction, Epping net | p.e. ± 0.86″ to ± 1.37″ | — |
| Single observed angle | — | ± 0.6″ to ± 1.0″ (1σ) in every survey except Cassini 1739 (± 8.8″) |
| Triangle closure, Eastern Oblique Arc | mean 1.09″ (46 triangles, Epping net) | σ ± 1.40″ over 267 triangles |
| Helmert's printed per-angle corrections | rms ± 0.65″ | ± 0.59″ from his closures |

The angles are the weak link: pointing noise (≈ ± 1.3″ per closure) is the only error close to the signal. Length and area errors are 10³–10⁵ times smaller. Large mountain-top triangles (up to 14 527 km², 75″ of excess) carry most of the result.

## Source check

45 triangles (135 angles) on five pages of Schott 1900 and 1902 were compared with the KML entry by entry: every value matches, including Schott's convention of carrying seconds above 60. Three triangles that do not close are misprints in Schott's own tables (one degree figure each), corrected from the printed side lengths by the law of cosines and kept. Two more (Helmert, Lambton) cannot be checked yet and are dropped.

## Contents

| Folder | What's in it |
| --- | --- |
| `index.html` | The report |
| `RULE.md` | The selection and weighting rule, written before the main run |
| `figures/` | Every chart in the report |
| `images/` | Plates from Schott 1900 and 1902: theodolites, heliotropes, mountain stations, towers, the Epping base |
| `supplemental/` | `Curvature_from_triangles.xlsx`: all 1 005 triangles; the main result, per-survey results and chain sums as live formulas (change the SNR cut or tier weights on `Inputs` and it recalculates); the source check with the law-of-cosines corrections |
| `sources/` | Page extracts from the scans. `*_annotated.pdf` carry the checks as PDF notes: green = verified, red = misprint |
| `data/` | The input KML, per-triangle results (`triangles_rule.csv`, `triangles_excess.csv`), fit and chain results as JSON |
| `scripts/` | Python that reproduces every number and chart; see `scripts/README.md` |
| [`model/`](https://github.com/rooster-ninja/TheoremaEgregium_CurvatureMeasuredFromSurveyorsTriangles/tree/main/model) | 3D model: North America on a dome built with this report's measured curvature (K = 2.450 × 10⁻¹⁴ m⁻²), with the US triangulation chains drawn on it, an ocean layer and the sky as a star plane. [Model page with rotatable viewer](https://rooster-ninja.github.io/TheoremaEgregium_CurvatureMeasuredFromSurveyorsTriangles/model/), [STL](https://github.com/rooster-ninja/TheoremaEgregium_CurvatureMeasuredFromSurveyorsTriangles/tree/main/model/NorthAmerica_EarthDome_K2.450.stl), SolidCraft file, an angle-deficit check of K, and where the population lives on the dome |
| [`theorem/`](https://github.com/rooster-ninja/TheoremaEgregium_CurvatureMeasuredFromSurveyorsTriangles/tree/main/theorem) | Background: the Theorema Egregium itself. History from Girard to Riemann, the name, what makes a theorem, sphere versus plane, and Earth-scale spherical-excess and circle-deficit tables for a distance-only test. [Read the page](https://rooster-ninja.github.io/TheoremaEgregium_CurvatureMeasuredFromSurveyorsTriangles/theorem/) |

## Sources

- Schott, C. A., *Geodesy: The Transcontinental Triangulation and the American Arc of the Parallel*, U.S. Coast and Geodetic Survey Special Publication No. 4, 1900.
- Schott, C. A., *Geodesy: The Eastern Oblique Arc of the United States and Osculating Spheroid*, U.S. Coast and Geodetic Survey Special Publication No. 7, 1902.
- Maclear, T. (ed. Airy), *Verification and Extension of La Caille's Arc of Meridian at the Cape of Good Hope*, 2 vols, 1866 (not yet checked against the KML).
- Helmert 1893; Lambton 1816–20; Roy 1790; Morris & Gill 1908; Cassini de Thury 1744; Brousseaud 1839, as cited in the KML.
- Gauss, *Disquisitiones generales circa superficies curvas* (1827), Arts. 12 and 20.
