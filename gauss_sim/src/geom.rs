//! Constant-curvature geometry from surface distances only.
//!
//! Every length is a geodesic length measured on the surface. K is the Gaussian
//! curvature. No angle and no radius is an input anywhere in this module.

use std::f64::consts::PI;

/// Σₖ (−z)ᵏ / (2k + first)!  (Stumpff series; `first` = 2 gives c₂, 3 gives c₃).
fn stumpff_series(z: f64, first: u32) -> f64 {
    let mut term = 1.0 / (1..=first).map(f64::from).product::<f64>();
    let mut sum = 0.0;
    for k in 0..14u32 {
        sum += term;
        let j = f64::from(2 * k + first);
        term *= -z / ((j + 1.0) * (j + 2.0));
    }
    sum
}

/// Stumpff c₂(z) = (1 − cos√z) / z, continued to z ≤ 0.
/// Half-angle form 2 sin²(√z/2)/z avoids the cancellation in 1 − cos√z.
pub fn c2(z: f64) -> f64 {
    if z.abs() < 1.0 {
        stumpff_series(z, 2)
    } else if z > 0.0 {
        2.0 * (z.sqrt() / 2.0).sin().powi(2) / z
    } else {
        2.0 * ((-z).sqrt() / 2.0).sinh().powi(2) / -z
    }
}

/// Stumpff c₃(z) = (√z − sin√z) / z^{3/2}, continued to z ≤ 0.
pub fn c3(z: f64) -> f64 {
    if z.abs() < 1.0 {
        stumpff_series(z, 3)
    } else if z > 0.0 {
        let s = z.sqrt();
        (s - s.sin()) / (z * s)
    } else {
        let s = (-z).sqrt();
        (s.sinh() - s) / (-z * s)
    }
}

/// Generalized sine snₖ(x): sin(√K x)/√K, x, sinh(√−K x)/√−K.
pub fn sn(k: f64, x: f64) -> f64 {
    x * (1.0 - k * x * x * c3(k * x * x))
}

/// Generalized cosine csₖ(x): cos(√K x), 1, cosh(√−K x).
pub fn cs(k: f64, x: f64) -> f64 {
    1.0 - k * x * x * c2(k * x * x)
}

/// Helper for the area formula: 2 tan(√K y/2)/√K, y, 2 tanh(√−K y/2)/√−K.
pub fn tn(k: f64, y: f64) -> f64 {
    if k > 0.0 {
        let s = k.sqrt();
        2.0 * (s * y / 2.0).tan() / s
    } else if k < 0.0 {
        let s = (-k).sqrt();
        2.0 * (s * y / 2.0).tanh() / s
    } else {
        y
    }
}

/// Most negative K·L² searched. Beyond it Gₖ is a difference of e^{2√|K|L}-sized terms
/// and loses its sign in double precision.
const K_FLOOR: f64 = -400.0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CutError {
    /// The three sides fail the triangle inequality.
    NotATriangle,
    /// The foot is not strictly inside the base (0 < m < ℓ).
    FootOffBase,
    /// d is too short for every K ≥ K_FLOOR/L².
    TooShort,
    /// d is too long for every K that holds the triangle.
    TooLong,
}

impl std::fmt::Display for CutError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str(match self {
            CutError::NotATriangle => "sides fail the triangle inequality",
            CutError::FootOffBase => "cut end off its side (need 0 < m < ℓ; cross: 0 < x ≤ u, 0 < y ≤ v)",
            CutError::TooShort => "d too short: no constant-K surface with R > L/20 fits",
            CutError::TooLong => "d too long: no constant-K surface fits",
        })
    }
}

/// One solver evaluation: the bracket [lo, hi], the K tried and Gₖ there.
#[derive(Clone, Copy, Debug)]
pub struct Step {
    pub stage: &'static str,
    pub lo: f64,
    pub hi: f64,
    pub k: f64,
    pub g: f64,
}

/// One cut in triangle PQQ′, with `base` = ℓ = |QQ′|, `near` = u = |PQ|, `far` = v = |PQ′|.
///
/// Foot cut (`y` = None): apex P to foot D on the base; `m` = |QD| along the base from Q,
/// `d` = |PD|, n = ℓ − m.
/// Cross cut (`y` = Some): X on PQ with `m` = x = |PX|, Y on PQ′ with y = |PY|, `d` = |XY|.
#[derive(Clone, Copy, Debug)]
pub struct Cut {
    pub base: f64,
    pub near: f64,
    pub far: f64,
    pub m: f64,
    pub y: Option<f64>,
    pub d: f64,
}

/// cos θ of the angle between sides s₁ and s₂ with opposite side o, from the law of cosines (1):
/// [csₖ(o) − csₖ(s₁)csₖ(s₂)] / (K snₖ(s₁)snₖ(s₂)), in a form exact at K = 0.
pub fn cos_angle(k: f64, s1: f64, s2: f64, o: f64) -> f64 {
    let q = |x: f64| x * x * c2(k * x * x);
    (q(s1) + q(s2) - q(o) - k * q(s1) * q(s2)) / (sn(k, s1) * sn(k, s2))
}

impl Cut {
    pub fn n(&self) -> f64 {
        self.base - self.m
    }

    /// Largest length in the cut; the natural length scale L.
    pub fn scale(&self) -> f64 {
        self.base.max(self.near).max(self.far).max(self.d)
    }

    /// Position of the cut along its side: t = m/ℓ (foot) or x/u (cross).
    pub fn t(&self) -> f64 {
        self.m / if self.y.is_some() { self.near } else { self.base }
    }

    pub fn check(&self) -> Result<(), CutError> {
        let (l, u, v) = (self.base, self.near, self.far);
        // Squares of the lengths must stay finite.
        if ![l, u, v, self.m, self.d, self.y.unwrap_or(1.0)].iter().all(|x| x.is_finite() && x.abs() < 1e150) {
            return Err(CutError::NotATriangle);
        }
        if !(l > 0.0 && u > 0.0 && v > 0.0 && l < u + v && u < l + v && v < l + u) {
            return Err(CutError::NotATriangle);
        }
        let on_sides = match self.y {
            None => self.m > 0.0 && self.m < l,
            Some(y) => self.m > 0.0 && self.m <= u && y > 0.0 && y <= v && !(self.m == u && y == v),
        };
        if !on_sides {
            return Err(CutError::FootOffBase);
        }
        if !(self.d > 0.0) {
            return Err(CutError::TooShort);
        }
        Ok(())
    }

    /// Cut function Gₖ = [csₖ(v) snₖ(m) + csₖ(u) snₖ(n) − csₖ(d) snₖ(ℓ)] / K  (equation 5).
    /// Written with Stumpff functions so it stays exact through K = 0, where Gₖ = Δ/2.
    /// Cross cut: Gₖ = cos θ(u, v; ℓ) − cos θ(x, y; d), the angle at P eliminated (equation 2b).
    pub fn g(&self, k: f64) -> f64 {
        let (l, u, v, m, d) = (self.base, self.near, self.far, self.m, self.d);
        if let Some(y) = self.y {
            return cos_angle(k, u, v, l) - cos_angle(k, m, y, d);
        }
        let n = self.n();
        let c3l = |x: f64| x * x * x * c3(k * x * x);
        let c2l = |x: f64| x * x * c2(k * x * x);
        c3l(l) - c3l(m) - c3l(n) - c2l(v) * sn(k, m) - c2l(u) * sn(k, n) + c2l(d) * sn(k, l)
    }

    /// Flat-plane prediction of the cut length.
    /// Foot (Stewart 1746): d₀² = (v²m + u²n)/ℓ − mn. Cross: d₀² = x² + y² − xy(u² + v² − ℓ²)/(uv).
    pub fn d_flat(&self) -> f64 {
        let (l, u, v, m) = (self.base, self.near, self.far, self.m);
        match self.y {
            None => ((v * v * m + u * u * self.n()) / l - m * self.n()).sqrt(),
            Some(y) => (m * m + y * y - m * y * (u * u + v * v - l * l) / (u * v)).max(0.0).sqrt(),
        }
    }

    /// Defect Δ = d² − d₀². Zero on a plane.
    pub fn defect(&self) -> f64 {
        let (l, u, v, m, n, d) = (self.base, self.near, self.far, self.m, self.n(), self.d);
        match self.y {
            None => (l * (d * d + m * n) - v * v * m - u * u * n) / l,
            Some(_) => d * d - self.d_flat().powi(2),
        }
    }

    /// Ratio ρ = d / d₀. ρ = 1 flat, ρ > 1 for K > 0, ρ < 1 for K < 0.
    pub fn rho(&self) -> f64 {
        self.d / self.d_flat()
    }

    /// Largest K that still holds the triangle and the cut: perimeter·√K < 2π, d·√K < π.
    pub fn k_max(&self) -> f64 {
        let half_p = (self.base + self.near + self.far) / 2.0;
        (PI / half_p.max(self.d)).powi(2)
    }

    /// Solve Gₖ = 0 for K. Gₖ > 0 below the root and < 0 above it, so bisection is exact.
    pub fn solve_k(&self) -> Result<f64, CutError> {
        self.solve_k_traced(&mut |_| ())
    }

    /// `solve_k`, reporting every evaluation as (stage, lo, hi, K tried, Gₖ there).
    pub fn solve_k_traced(&self, log: &mut dyn FnMut(Step)) -> Result<f64, CutError> {
        self.check()?;
        let g0 = self.g(0.0);
        log(Step { stage: "start", lo: 0.0, hi: 0.0, k: 0.0, g: g0 });
        if g0 == 0.0 {
            return Ok(0.0);
        }
        let l2 = self.scale().powi(2);
        let (mut lo, mut hi);
        if g0 > 0.0 {
            lo = 0.0;
            hi = self.k_max() * (1.0 - 1e-12);
            let gh = self.g(hi);
            log(Step { stage: "bracket", lo, hi, k: hi, g: gh });
            if gh > 0.0 {
                return Err(CutError::TooLong);
            }
        } else {
            hi = 0.0;
            lo = -1.0 / l2;
            loop {
                let gl = self.g(lo);
                log(Step { stage: "bracket", lo, hi, k: lo, g: gl });
                if gl >= 0.0 {
                    break;
                }
                if lo * l2 <= K_FLOOR {
                    return Err(CutError::TooShort);
                }
                lo = (2.0 * lo).max(K_FLOOR / l2);
            }
        }
        for _ in 0..400 {
            let mid = 0.5 * (lo + hi);
            if mid <= lo || mid >= hi {
                break;
            }
            let gm = self.g(mid);
            log(Step { stage: "bisect", lo, hi, k: mid, g: gm });
            if gm > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(0.5 * (lo + hi))
    }

    /// Forward problem: the cut length d on a surface of curvature K (d is ignored on input).
    pub fn length_on(&self, k: f64) -> f64 {
        let mut c = *self;
        let (mut lo, mut hi) = match self.y {
            None => (0.0, self.m + self.near),
            Some(y) => ((self.m - y).abs(), self.m + y),
        };
        if k > 0.0 {
            hi = hi.min(PI / k.sqrt());
        }
        for _ in 0..200 {
            c.d = 0.5 * (lo + hi);
            // Gₖ grows with d.
            if c.g(k) > 0.0 {
                hi = c.d;
            } else {
                lo = c.d;
            }
        }
        0.5 * (lo + hi)
    }
}

/// Angle (radians) between sides s₁ and s₂ of a geodesic triangle with opposite side o.
/// Half-angle form of the law of cosines (1), with s the semi-perimeter:
/// sin²(θ/2) = snₖ(s − s₁)·snₖ(s − s₂) / (snₖ(s₁)·snₖ(s₂)). No cancellation for any K.
pub fn angle(k: f64, s1: f64, s2: f64, o: f64) -> f64 {
    let s = (s1 + s2 + o) / 2.0;
    let h = sn(k, s - s1) * sn(k, s - s2) / (sn(k, s1) * sn(k, s2));
    2.0 * h.clamp(0.0, 1.0).sqrt().asin()
}

/// Area 𝒜 of a geodesic triangle with sides a, b, c on curvature K
/// (L'Huilier form; Heron at K = 0). Angle excess ε = K𝒜 (Gauss–Bonnet).
pub fn area(k: f64, a: f64, b: f64, c: f64) -> f64 {
    let s = (a + b + c) / 2.0;
    let h = (tn(k, s) * tn(k, s - a) * tn(k, s - b) * tn(k, s - c)).max(0.0).sqrt();
    if k == 0.0 {
        h
    } else {
        4.0 * (k.abs() * h / 4.0).atan() / k.abs()
    }
}

/// A point in the model surface of curvature K, with vertex A at the pole/origin.
/// K > 0: sphere |x| = R in ℝ³. K < 0: hyperboloid x² + y² − z² = −R². K = 0: plane z = 0.
pub type P3 = [f64; 3];

#[derive(Clone, Copy, Debug)]
pub struct Model {
    pub k: f64,
}

impl Model {
    fn r(&self) -> f64 {
        1.0 / self.k.abs().sqrt()
    }

    /// Point at geodesic distance `dist` from the pole, in direction `az` (radians).
    pub fn polar_point(&self, dist: f64, az: f64) -> P3 {
        let (sa, ca) = az.sin_cos();
        if self.k > 0.0 {
            let r = self.r();
            let t = dist / r;
            [r * t.sin() * ca, r * t.sin() * sa, r * t.cos()]
        } else if self.k < 0.0 {
            let r = self.r();
            let t = dist / r;
            [r * t.sinh() * ca, r * t.sinh() * sa, r * t.cosh()]
        } else {
            [dist * ca, dist * sa, 0.0]
        }
    }

    /// Geodesic distance from the pole and azimuth (azimuthal-equidistant map about A).
    pub fn polar(&self, p: P3) -> (f64, f64) {
        let rho = p[0].hypot(p[1]);
        let az = p[1].atan2(p[0]);
        let dist = if self.k > 0.0 {
            self.r() * rho.atan2(p[2])
        } else if self.k < 0.0 {
            self.r() * (rho / self.r()).asinh()
        } else {
            rho
        };
        (dist, az)
    }

    /// Point a fraction t of the way along the geodesic from p to q.
    pub fn lerp(&self, p: P3, q: P3, t: f64) -> P3 {
        if self.k == 0.0 {
            return std::array::from_fn(|i| p[i] + t * (q[i] - p[i]));
        }
        let r = self.r();
        let (u, v): (P3, P3) = (p.map(|x| x / r), q.map(|x| x / r));
        let w: P3 = std::array::from_fn(|i| u[i] - v[i]);
        // Chord length gives θ without cancellation: 2 sin(θ/2) or 2 sinh(θ/2).
        let (theta, s): (f64, fn(f64) -> f64) = if self.k > 0.0 {
            (2.0 * ((w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt() / 2.0).min(1.0).asin(), f64::sin)
        } else {
            let lw = (w[0] * w[0] + w[1] * w[1] - w[2] * w[2]).max(0.0).sqrt();
            // The chord form cancels for long arcs; there cosh θ = −⟨u, v⟩ is exact enough.
            let theta = if lw < 1.0 { 2.0 * (lw / 2.0).asinh() } else { (u[2] * v[2] - u[0] * v[0] - u[1] * v[1]).max(1.0).acosh() };
            (theta, f64::sinh)
        };
        if theta < 1e-15 {
            return p;
        }
        let (fa, fb) = (s((1.0 - t) * theta) / s(theta), s(t * theta) / s(theta));
        std::array::from_fn(|i| r * (fa * u[i] + fb * v[i]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rel(a: f64, b: f64) -> f64 {
        (a - b).abs() / b.abs().max(1e-300)
    }

    #[test]
    fn stumpff_branches_join() {
        for z in [0.9999999f64, -0.9999999, 1.0000001, -1.0000001] {
            let direct2 = if z > 0.0 { (1.0 - z.sqrt().cos()) / z } else { ((-z).sqrt().cosh() - 1.0) / -z };
            assert!(rel(stumpff_series(z, 2), direct2) < 1e-13);
            assert!(rel(c2(z), direct2) < 1e-13);
        }
        assert_eq!(c2(0.0), 0.5);
        assert!(rel(c3(0.0), 1.0 / 6.0) < 1e-16);
    }

    #[test]
    fn g_at_zero_is_half_defect() {
        let cut = Cut { base: 3.0, near: 2.5, far: 2.0, m: 1.2, y: None, d: 1.9 };
        assert!(rel(cut.g(0.0), cut.base * cut.defect() / 2.0) < 1e-14);
    }

    /// Exact octant of the unit sphere: all sides π/2, median from A has length π/2.
    #[test]
    fn octant_recovers_k_one() {
        let h = PI / 2.0;
        let cut = Cut { base: h, near: h, far: h, m: h / 2.0, y: None, d: h };
        assert!(cut.rho() > 1.0);
        let k = cut.solve_k().unwrap();
        assert!(rel(k, 1.0) < 1e-12, "k = {k}");
    }

    #[test]
    fn plane_gives_zero() {
        // 3-4-5 right triangle, apex at the right angle; cut to the hypotenuse midpoint has d = 2.5.
        let cut = Cut { base: 5.0, near: 3.0, far: 4.0, m: 2.5, y: None, d: 2.5 };
        assert!(cut.defect().abs() < 1e-12);
        assert!(cut.solve_k().unwrap().abs() < 1e-15);
    }

    #[test]
    fn forward_then_inverse_round_trip() {
        for &k in &[2.0, 0.3, 1e-4, 1e-9, 0.0, -1e-9, -1e-4, -0.3, -2.0] {
            for &(l, u, v, t) in &[(1.0, 0.8, 0.7, 0.5), (0.9, 1.1, 0.5, 0.2), (1.2, 0.6, 0.9, 0.85)] {
                let mut cut = Cut { base: l, near: u, far: v, m: t * l, y: None, d: 0.0 };
                cut.d = cut.length_on(k);
                let back = cut.solve_k().unwrap();
                // Double-precision floor: |δK|·L² ≈ 1e-14 near K = 0.
                assert!((back - k).abs() < 1e-9 * k.abs() + 1e-13, "k={k} back={back}");
                if k != 0.0 {
                    assert_eq!(cut.rho() > 1.0, k > 0.0);
                }
            }
        }
    }

    /// Distances measured directly in the embedded model agree with the cut equation.
    #[test]
    fn model_embedding_matches_cut_equation() {
        // −60 makes |BC| ≈ 7 R: the long-arc branch of lerp.
        for &k in &[1.7, -1.7, -60.0] {
            let model = Model { k };
            let (a, b, c) = (1.0, 0.9, 0.7);
            let alpha = angle(k, b, c, a);
            let pb = model.polar_point(c, 0.0);
            let pc = model.polar_point(b, alpha);
            let pd = model.lerp(pb, pc, 0.3);
            // distance |AD| read off the pole map equals the forward cut length.
            let (d, _) = model.polar(pd);
            let cut = Cut { base: a, near: c, far: b, m: 0.3 * a, y: None, d: 0.0 };
            assert!(rel(d, cut.length_on(k)) < 1e-9, "k={k}: {d} vs {}", cut.length_on(k));
            assert!(rel(model.polar(pc).0, b) < 1e-12);
        }
    }

    #[test]
    fn gauss_bonnet_angle_excess() {
        for &k in &[1.3, 0.01, -0.01, -1.3, -50.0, -300.0] {
            let (a, b, c) = (1.0, 0.8, 0.9);
            let sum = angle(k, b, c, a) + angle(k, c, a, b) + angle(k, a, b, c);
            let excess = k * area(k, a, b, c);
            assert!(((sum - PI) - excess).abs() < 1e-12, "k={k}: {} vs {excess}", sum - PI);
        }
    }

    /// Strong negative curvature down to the search floor is still recovered.
    #[test]
    fn saddle_near_floor() {
        for &k in &[-50.0, -350.0] {
            let mut cut = Cut { base: 1.0, near: 0.9, far: 0.8, m: 0.5, y: None, d: 0.0 };
            cut.d = cut.length_on(k);
            assert!(rel(cut.solve_k().unwrap(), k) < 1e-6);
        }
        let cut = Cut { base: 1.0, near: 0.9, far: 0.8, m: 0.5, y: None, d: 0.39 };
        assert_eq!(cut.solve_k(), Err(CutError::TooShort));
    }

    /// Cross cut joining the midpoints of the two legs from P (midline): d₀ = ℓ/2 on a plane,
    /// and the sphere and saddle values come back from the solver.
    #[test]
    fn cross_cut_round_trip() {
        let mut cut = Cut { base: 0.55, near: 1.0, far: 1.0, m: 0.5, y: Some(0.5), d: 0.0 };
        assert!(rel(cut.d_flat(), 0.275) < 1e-14);
        for &k in &[2.0, 1e-6, 0.0, -1e-6, -2.0] {
            cut.d = cut.length_on(k);
            let back = cut.solve_k().unwrap();
            assert!((back - k).abs() < 1e-8 * k.abs() + 1e-12, "k={k} back={back}");
            if k != 0.0 {
                assert_eq!(cut.rho() > 1.0, k > 0.0);
            }
        }
    }

    /// Sphere midline check against the law of cosines on a sphere of radius 600.
    #[test]
    fn cross_cut_on_sphere() {
        let (r, p, th): (f64, f64, f64) = (600.0, 155.0, 0.56);
        let leg = |s: f64| r * ((s / r).cos().powi(2) + (s / r).sin().powi(2) * th.cos()).acos();
        let cut = Cut { base: leg(p), near: p, far: p, m: p / 2.0, y: Some(p / 2.0), d: leg(p / 2.0) };
        assert!(rel(cut.solve_k().unwrap(), r.powi(-2)) < 1e-6);
    }
}
