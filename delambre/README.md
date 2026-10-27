# Delambre and Méchain: the meridian arc from Dunkerque to Barcelona

### ▶ [Read the Évaux to Carcassonne report](https://alanspaceaudits.github.io/theorema-egregium/delambre/)

Two determinations from the printed record of the French meridian arc (1792 to 1798).

## 1. Évaux to Carcassonne: the sky against the ground

**Claim checked.** "From Évaux to Carcassonne the physically measured distance is 350 km. Astronomical latitude says it is 330 km."

**Result.** The printed ground arc is 168,846.7 toises = 329.088 km, not 350 km. The sky and the ground agree to within 117 to 250 m (3.8″ to 8.1″ of latitude), depending on which degree length is used. A 20 km mismatch needs an error of 648″ (10.8′) in the latitudes. Each end latitude is the mean of hundreds of star observations with the Borda repeating circle.

The full report with the cited pages is `index.html`.

## 2. Curvature from the observed angle sums

The fit uses the 95 main triangles that have all three angles observed. Each area is Heron's formula on the three printed sides. No radius is assumed.

| Fit | K | R |
|---|---|---|
| Through zero (observed sums) | 1.319e-13 per square toise = 3.47e-14 per square metre | 2,753,661 toises = 5,367 km (1σ 5,166 to 5,594 km) |
| With a constant bias per triangle (bias +0.24″) | | 5,634 km |
| Delambre's printed excess, for comparison | | 6,383 km |

The through-zero fit is 12.6σ from K = 0 (a plane). The scatter of the sums about the fit is 1.76″. The full output is in `analysis/curvature_fit_output.txt`.

## Files

| Path | Contents |
|---|---|
| `index.html` | Évaux to Carcassonne report. One self-contained file. |
| `analysis/review.json` | The 115 main triangles: observed angles, excess, spherical angles, register pages, sums and sides. 93 station latitudes from vol. 3, Table V. |
| `analysis/astro.json` | Latitude and azimuth stations, star examples, partial arcs and degree lengths. |
| `analysis/observations.json` | The vol. 1 field registers for 96 stations: 329 angle entries and 217 sighted lines. |
| `analysis/work/t1_transcription.txt` | Transcription of the vol. 1 Tableau général des triangles. |
| `analysis/check.py` | Checks of the transcription. Output: `check_output.txt`. |
| `analysis/curvature_fit.py` | The curvature fit. Output: `curvature_fit_output.txt`. |

To run both scripts, go to `analysis/` and run `python3 check.py` and `python3 curvature_fit.py`. They need only the Python standard library.

## What the printed angles are

Vol. 1, pp. 513 to 542, prints the Tableau général des triangles. Each triangle angle appears in five columns:

1. **Angles observés.** The commission's value from the register series. It is reduced to the centre and to the horizon. No closing correction is applied.
2. **Excès sphérique.** One value per angle, from tables built on a degree of 57,020 toises (p. 147).
3. **Angles sphériques.** The observed angle plus one third of the closing error.
4. **Angles des cordes.**
5. **Angles moyens.**

Vol. 2, pp. 800 to 839, prints the Tableau complet with the sides. Its angles are the vol. 1 spherical angles. In triangles 44 to 100, Delambre changed them by 0.01″ to 0.1″ so that the Melun and Perpignan bases agree (pp. 704 to 705).

The observed sums in this folder come from the vol. 1 observed column only.

Original (French), vol. 1, p. 511:
> La quatrième offre les angles tels qu'ils ont été arrêtés par la commission spéciale chargée d'examiner nos registres.

Translation (English):
> The fourth gives the angles as the special commission charged with examining our registers fixed them.

Original (French), vol. 2, p. 700:
> les angles sphériques tels qu'on les a donnés tom. I, p. 513 et suiv.; c'est-à-dire corrigés du tiers de la somme des trois erreurs.

Translation (English):
> the spherical angles as given in vol. I, p. 513 and following; that is, corrected by one third of the sum of the three errors.

## Transcription

- **Triangles.** All 115 main triangles, Dunkerque to Montjoui. 95 have all three angles observed. The other 20 have one angle concluded ('*'), so their sum is forced and they are not in the fit.
- **Closure.** Mean closure (observed sum minus 180° minus the printed excess) is +0.64″, standard deviation 1.76″.
- **Sides.** 231 side entries in toises from vol. 2. Metres = toises × 1.949036591. This factor is 10^7 / 5,130,740, from the quarter meridian that the commission adopted (vol. 3, p. 107).
- **Bases.** Melun (Lieursaint to Melun, 6,075.90 toises) and Perpignan (Vernet to Salces, 6,006.25 toises), vol. 2, p. 56.
- **Stations.** Latitudes for 93 stations from vol. 3, Table V, pp. 266 to 268. Table V prints longitudes only from the meridian of Dunkerque, so `lon` is null.

## Checks (check.py)

1. **Sums.** The recomputed observed sums match the printed sum line in 114 of 115 triangles. In triangle 60 the printed angles add to 179° 0′ 0.0″ but the sum line prints 180° 0′ 0.0″. That triangle has a concluded angle and is not in the fit.
2. **Sine rule.** The spread of side / sin(angle − E/3) is within 1e-5 in 110 of 115 triangles. In triangles 29, 44, 61 and 87 the printed arc disagrees with the printed log side, and the log side fits the angles. Triangle 68 prints its sides to 0.1 toise.
3. **Shared sides.** 13 sides are printed in two triangles. They differ by 0.02 to 0.81 toise, from different routes through the net.

## Gaps

- No Greenwich or Paris longitudes in Table V.
- The secondary triangles are not transcribed: vol. 1, pp. 542 to 551, including the Barcelona end, and vol. 2 from p. 840.
- No probable error of an angle is printed.

## Sources

The full source list, with every cited page, is in [`SOURCES.md`](SOURCES.md). The three volumes are the ETH-Bibliothek Zürich copy on e-rara: https://doi.org/10.3931/e-rara-11803

- Delambre, J.-B. J., & Méchain, P. F. A. (1806). *Base du système métrique décimal, ou mesure de l'arc du méridien compris entre les parallèles de Dunkerque et Barcelone, exécutée en 1792 et années suivantes* (Vol. 1). Paris: Baudouin.
- Delambre, J.-B. J., & Méchain, P. F. A. (1807). *Base du système métrique décimal* (Vol. 2). Paris: Baudouin.
- Delambre, J.-B. J., & Méchain, P. F. A. (1810). *Base du système métrique décimal* (Vol. 3). Paris: Baudouin.
