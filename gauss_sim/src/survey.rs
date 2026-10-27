//! Measurement sets, curvature analysis across cuts, and the simulator.

use crate::geom::{angle, area, Cut, CutError};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Vertex {
    A,
    B,
    C,
}

impl Vertex {
    pub const ALL: [Vertex; 3] = [Vertex::A, Vertex::B, Vertex::C];

    pub fn name(self) -> &'static str {
        match self {
            Vertex::A => "A",
            Vertex::B => "B",
            Vertex::C => "C",
        }
    }

    /// Base of a cut from this apex, written from Q to Q′ (m is measured from Q).
    pub fn base_name(self) -> &'static str {
        match self {
            Vertex::A => "BC",
            Vertex::B => "CA",
            Vertex::C => "AB",
        }
    }

    pub fn index(self) -> usize {
        self as usize
    }
}

/// Measured sides: a = |BC|, b = |CA|, c = |AB|.
#[derive(Clone, Copy, Debug)]
pub struct Triangle {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

impl Triangle {
    /// (ℓ, u, v) for a cut from `apex` P, with (P, Q, Q′) cyclic in (A, B, C):
    /// ℓ = |QQ′|, u = |PQ|, v = |PQ′|.
    pub fn frame(&self, apex: Vertex) -> (f64, f64, f64) {
        match apex {
            Vertex::A => (self.a, self.c, self.b),
            Vertex::B => (self.b, self.a, self.c),
            Vertex::C => (self.c, self.b, self.a),
        }
    }

    pub fn cut(&self, c: &CutInput) -> Cut {
        let (base, near, far) = self.frame(c.apex);
        Cut { base, near, far, m: c.m, y: c.y, d: c.d }
    }

    pub fn scale(&self) -> f64 {
        self.a.max(self.b).max(self.c)
    }
}

/// One measured cut from apex P.
/// Foot cut (`y` = None): m = |QD| along the base, d = |PD|.
/// Cross cut (`y` = Some): m = x = |PX| on PQ, y = |PY| on PQ′, d = |XY|.
#[derive(Clone, Copy, Debug)]
pub struct CutInput {
    pub apex: Vertex,
    pub m: f64,
    pub y: Option<f64>,
    pub d: f64,
}

impl CutInput {
    /// The measured lengths of this cut, in a fixed order (m, [y,] d).
    fn lengths(&self) -> Vec<f64> {
        let mut v = vec![self.m];
        v.extend(self.y);
        v.push(self.d);
        v
    }

    fn with_lengths(&self, x: &[f64]) -> CutInput {
        let y = self.y.map(|_| x[1]);
        CutInput { apex: self.apex, m: x[0], y, d: x[x.len() - 1] }
    }
}

pub struct CutResult {
    pub cut: Cut,
    pub apex: Vertex,
    pub k: Result<f64, CutError>,
    /// σ(Kⱼ) from σ on the cut's lengths (ℓ, u, v, m, [y,] d).
    pub sigma_k: f64,
    pub sigma_rho: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Verdict {
    Flat,
    Positive,
    Negative,
}

pub struct Analysis {
    pub cuts: Vec<CutResult>,
    /// K*: combined curvature, inverse-variance mean of the Kⱼ.
    pub k_star: Option<f64>,
    pub sigma_k_star: Option<f64>,
    /// χ²/dof of the Kⱼ about K* (needs σ > 0 and ≥ 2 cuts).
    pub chi2_dof: Option<f64>,
    pub verdict: Option<Verdict>,
    /// Interior angles α, β, γ at A, B, C (radians), derived on K*.
    pub angles: Option<[f64; 3]>,
    pub area: Option<f64>,
    pub triangle_ok: bool,
}

/// σ of f(x) when every entry of x carries an independent error σ (central differences,
/// one-sided where a shift leaves the valid domain). NaN if f fails on both sides.
fn propagate(x: &[f64], sigma: f64, f: impl Fn(&[f64]) -> Option<f64>) -> f64 {
    if sigma <= 0.0 {
        return 0.0;
    }
    let h = 1e-6 * x.iter().fold(0.0, |m: f64, v| m.max(v.abs()));
    let Some(f0) = f(x) else { return f64::NAN };
    let mut var = 0.0;
    for i in 0..x.len() {
        let shift = |s: f64| {
            let mut y = x.to_vec();
            y[i] += s;
            f(&y)
        };
        let slope = match (shift(h), shift(-h)) {
            (Some(p), Some(q)) => (p - q) / (2.0 * h),
            (Some(p), None) => (p - f0) / h,
            (None, Some(q)) => (f0 - q) / h,
            (None, None) => return f64::NAN,
        };
        var += (slope * sigma).powi(2);
    }
    var.sqrt()
}

fn cut_vec(c: &Cut) -> Vec<f64> {
    let mut v = vec![c.base, c.near, c.far, c.m];
    v.extend(c.y);
    v.push(c.d);
    v
}

fn vec_cut(x: &[f64]) -> Cut {
    Cut { base: x[0], near: x[1], far: x[2], m: x[3], y: (x.len() == 6).then(|| x[4]), d: x[x.len() - 1] }
}

pub fn analyze(tri: &Triangle, inputs: &[CutInput], sigma: f64) -> Analysis {
    let probe = Cut { base: tri.a, near: tri.c, far: tri.b, m: tri.a / 2.0, y: None, d: tri.a };
    let triangle_ok = probe.check().is_ok();

    let cuts: Vec<CutResult> = inputs
        .iter()
        .map(|ci| {
            let cut = tri.cut(ci);
            let k = cut.solve_k();
            let (sigma_k, sigma_rho) = if k.is_ok() {
                let x = cut_vec(&cut);
                (
                    propagate(&x, sigma, |y| vec_cut(y).solve_k().ok()),
                    propagate(&x, sigma, |y| Some(vec_cut(y).rho())),
                )
            } else {
                (0.0, 0.0)
            };
            CutResult { cut, apex: ci.apex, k, sigma_k, sigma_rho }
        })
        .collect();

    // Cuts that solved: (index into inputs, Kⱼ, weight).
    let weighted = sigma > 0.0 && cuts.iter().all(|c| c.k.is_err() || (c.sigma_k.is_finite() && c.sigma_k > 0.0));
    let good: Vec<(usize, f64, f64)> = cuts
        .iter()
        .enumerate()
        .filter_map(|(i, c)| c.k.ok().map(|k| (i, k, if weighted { c.sigma_k.powi(-2) } else { 1.0 })))
        .collect();
    let nk = good.len();
    let wsum: f64 = good.iter().map(|g| g.2).sum();
    let k_star = (nk > 0).then(|| good.iter().map(|g| g.1 * g.2).sum::<f64>() / wsum);

    // σ(K*): K* as one function of every measured length (a, b, c, then each solved cut's lengths),
    // with the weights held fixed. Shared sides enter once, so their correlation is counted.
    let sigma_k_star = k_star.and_then(|k| {
        if sigma > 0.0 {
            let mut x = vec![tri.a, tri.b, tri.c];
            x.extend(good.iter().flat_map(|g| inputs[g.0].lengths()));
            let f = |y: &[f64]| {
                let t = Triangle { a: y[0], b: y[1], c: y[2] };
                let (mut acc, mut at) = (0.0, 3);
                for g in &good {
                    let ci = &inputs[g.0];
                    let n = ci.lengths().len();
                    acc += g.2 * t.cut(&ci.with_lengths(&y[at..at + n])).solve_k().ok()?;
                    at += n;
                }
                Some(acc / wsum)
            };
            Some(propagate(&x, sigma, f)).filter(|s| s.is_finite())
        } else {
            // Exact data: the scatter of the Kⱼ is the only spread available.
            (nk > 1).then(|| {
                let var = good.iter().map(|g| (g.1 - k).powi(2)).sum::<f64>() / (nk - 1) as f64;
                (var / nk as f64).sqrt()
            })
        }
    });
    let chi2_dof = k_star.filter(|_| weighted && nk > 1).map(|k| {
        good.iter().map(|g| (g.1 - k).powi(2) * g.2).sum::<f64>() / (nk - 1) as f64
    });

    let l2 = tri.scale().powi(2);
    let verdict = k_star.map(|k| {
        // Floating-point noise (|K|·L² ≈ 1e-13) never counts as curvature.
        let tol = sigma_k_star.map_or(0.0, |s| 2.0 * s).max(1e-12 / l2);
        if k.abs() <= tol {
            Verdict::Flat
        } else if k > 0.0 {
            Verdict::Positive
        } else {
            Verdict::Negative
        }
    });

    let geo = k_star.filter(|&k| triangle_ok && (k <= 0.0 || (tri.a + tri.b + tri.c) * k.sqrt() < 2.0 * std::f64::consts::PI));
    let angles = geo.map(|k| [angle(k, tri.b, tri.c, tri.a), angle(k, tri.c, tri.a, tri.b), angle(k, tri.a, tri.b, tri.c)]);
    let area = geo.map(|k| area(k, tri.a, tri.b, tri.c));

    Analysis { cuts, k_star, sigma_k_star, chi2_dof, verdict, angles, area, triangle_ok }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Surface {
    Plane,
    Sphere,
    Saddle,
}

/// Simulator: a triangle with true sides on a surface of known K.
/// Cuts are (apex, cross, t): foot at m = tℓ, or cross cut at x = tu, y = tv.
#[derive(Clone, Debug)]
pub struct Sim {
    pub surface: Surface,
    /// R of the sphere, or pseudo-radius 1/√(−K) of the saddle.
    pub radius: f64,
    pub tri: Triangle,
    pub cuts: Vec<(Vertex, bool, f64)>,
    /// σ: standard deviation of the noise on each length. 0 gives exact data.
    pub noise: f64,
    pub seed: u64,
}

impl Sim {
    pub fn k_true(&self) -> f64 {
        match self.surface {
            Surface::Plane => 0.0,
            Surface::Sphere => self.radius.powi(-2),
            Surface::Saddle => -self.radius.powi(-2),
        }
    }

    /// True lengths on the model surface, then Gaussian noise on every measured length.
    pub fn measure(&self) -> Result<(Triangle, Vec<CutInput>), String> {
        if self.surface != Surface::Plane && !(self.radius > 0.0 && self.radius.is_finite()) {
            return Err("R_sim must be > 0".into());
        }
        let k = self.k_true();
        let t = self.tri;
        let probe = Cut { base: t.a, near: t.c, far: t.b, m: t.a / 2.0, y: None, d: t.a };
        probe.check().map_err(|e| e.to_string())?;
        if k > 0.0 && (t.a + t.b + t.c) * k.sqrt() >= 2.0 * std::f64::consts::PI {
            return Err("triangle does not fit on this sphere (perimeter ≥ 2πR)".into());
        }
        let mut rng = Rng(self.seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        let mut noisy = |x: f64| x + self.noise * rng.gauss();
        let cuts = self
            .cuts
            .iter()
            .map(|&(apex, cross, frac)| {
                let (base, near, far) = t.frame(apex);
                let (m, y) = if cross { (frac * near, Some(frac * far)) } else { (frac * base, None) };
                let cut = Cut { base, near, far, m, y, d: 0.0 };
                CutInput { apex, m, y, d: cut.length_on(k) }
            })
            .collect::<Vec<_>>();
        let tri = Triangle { a: noisy(t.a), b: noisy(t.b), c: noisy(t.c) };
        let cuts = cuts.into_iter().map(|c| CutInput { apex: c.apex, m: noisy(c.m), y: c.y.map(&mut noisy), d: noisy(c.d) }).collect();
        Ok((tri, cuts))
    }
}

/// SplitMix64 with Box–Muller.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    }

    fn gauss(&mut self) -> f64 {
        let u = self.next().max(1e-300);
        (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * self.next()).cos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sim(surface: Surface, radius: f64, side: f64, noise: f64) -> Sim {
        Sim {
            surface,
            radius,
            tri: Triangle { a: side, b: side * 0.9, c: side * 1.1 },
            cuts: Vertex::ALL.iter().flat_map(|&v| [(false, 0.25), (false, 0.5), (false, 0.75), (true, 0.5)].map(|(x, t)| (v, x, t))).collect(),
            noise,
            seed: 7,
        }
    }

    #[test]
    fn exact_sphere_all_cuts_agree() {
        let s = sim(Surface::Sphere, 6371.0, 300.0, 0.0);
        let (tri, cuts) = s.measure().unwrap();
        let an = analyze(&tri, &cuts, 0.0);
        for c in &an.cuts {
            let k = c.k.unwrap();
            assert!((k - s.k_true()).abs() < 1e-6 * s.k_true(), "{k} vs {}", s.k_true());
        }
        assert_eq!(an.verdict, Some(Verdict::Positive));
    }

    #[test]
    fn exact_saddle_and_plane() {
        let s = sim(Surface::Saddle, 50.0, 30.0, 0.0);
        let (tri, cuts) = s.measure().unwrap();
        let an = analyze(&tri, &cuts, 0.0);
        assert!((an.k_star.unwrap() / s.k_true() - 1.0).abs() < 1e-9);
        assert_eq!(an.verdict, Some(Verdict::Negative));

        let s = sim(Surface::Plane, 1.0, 30.0, 0.0);
        let (tri, cuts) = s.measure().unwrap();
        assert_eq!(analyze(&tri, &cuts, 0.0).verdict, Some(Verdict::Flat));
    }

    /// Noisy data: K* lands within a few σ of the truth.
    #[test]
    fn noisy_sphere_within_sigma() {
        let s = sim(Surface::Sphere, 6371.0, 500.0, 0.001);
        let (tri, cuts) = s.measure().unwrap();
        let an = analyze(&tri, &cuts, s.noise);
        let (k, sk) = (an.k_star.unwrap(), an.sigma_k_star.unwrap());
        assert!((k - s.k_true()).abs() < 4.0 * sk, "{k} ± {sk} vs {}", s.k_true());
        assert_eq!(an.verdict, Some(Verdict::Positive));
    }

    /// Exact plane data never reads as curved, for any triangle shape.
    #[test]
    fn exact_plane_is_flat() {
        for i in 0..200 {
            let mut s = sim(Surface::Plane, 1.0, 3.0, 0.0);
            s.tri.b = 2.0 + i as f64 * 0.02;
            let (tri, cuts) = s.measure().unwrap();
            assert_eq!(analyze(&tri, &cuts, 0.0).verdict, Some(Verdict::Flat), "b = {}", s.tri.b);
        }
    }

    /// Noisy plane data: the 2σ verdict calls "curved" about 5% of the time.
    #[test]
    fn flat_false_alarm_rate() {
        let mut curved = 0;
        for seed in 0..400 {
            let mut s = sim(Surface::Plane, 1.0, 1.0, 1e-3);
            s.seed = seed;
            let (tri, cuts) = s.measure().unwrap();
            curved += (analyze(&tri, &cuts, s.noise).verdict != Some(Verdict::Flat)) as u32;
        }
        assert!((8..=36).contains(&curved), "{curved}/400");
    }

    #[test]
    fn no_cut_no_angles() {
        let an = analyze(&Triangle { a: 1.0, b: 1.0, c: 1.0 }, &[], 0.0);
        assert!(an.verdict.is_none() && an.angles.is_none());
        let bad = Sim { radius: 0.0, ..sim(Surface::Saddle, 1.0, 1.0, 0.0) };
        assert!(bad.measure().is_err());
    }

    /// Midline cross cut on measured data: O = A, legs p = 155.477388, s₁ = 85.2099027, s₂ at p/2 = 42.9374544.
    #[test]
    fn midline_cross_cut_example() {
        let (p, s1, s2) = (155.477388, 85.2099027, 42.9374544);
        let tri = Triangle { a: s1, b: p, c: p };
        let cuts = [CutInput { apex: Vertex::A, m: p / 2.0, y: Some(p / 2.0), d: s2 }];
        let an = analyze(&tri, &cuts, 0.0);
        let r = 1.0 / an.k_star.unwrap().sqrt();
        assert!((r - 600.6123214782).abs() < 1e-4, "R = {r}");
        assert_eq!(an.verdict, Some(Verdict::Positive));
        let text = crate::report::build(&tri, &cuts, 0.0, "m", &an, None);
        assert!(text.contains("600.61"));
        println!("{text}");
    }
}
