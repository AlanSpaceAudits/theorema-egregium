//! GUI: simulator, measurement entry, results, plots and surface views.

use crate::geom::{Cut, Model, P3};
use crate::survey::{analyze, Analysis, CutInput, Sim, Surface, Triangle, Verdict, Vertex};
use eframe::egui::{self, Color32, Pos2, Stroke, Vec2};
use egui_plot::{HLine, Legend, Line, LineStyle, Plot, PlotPoints, Points};
use std::f64::consts::PI;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Tab {
    Results,
    Ratio,
    CutFunctions,
    Map,
    Sphere,
    Notation,
    Report,
}

const TABS: [(Tab, &str); 7] = [
    (Tab::Results, "Results"),
    (Tab::Ratio, "Ratio test ρ"),
    (Tab::CutFunctions, "Cut functions Gₖ"),
    (Tab::Map, "Intrinsic map"),
    (Tab::Sphere, "Sphere view"),
    (Tab::Notation, "Notation & method"),
    (Tab::Report, "Report"),
];

fn apex_color(v: Vertex) -> Color32 {
    [Color32::from_rgb(235, 110, 90), Color32::from_rgb(110, 200, 120), Color32::from_rgb(100, 160, 240)][v.index()]
}

/// Six significant figures; scientific notation outside 10⁻³ … 10⁶.
fn num(x: f64) -> String {
    if !x.is_finite() {
        return "n/a".into();
    }
    if x == 0.0 {
        return "0".into();
    }
    let a = x.abs();
    if (1e-3..1e6).contains(&a) {
        let dec = (5 - a.log10().floor() as i32).clamp(0, 9) as usize;
        format!("{x:.dec$}")
    } else {
        format!("{x:.5e}")
    }
}

struct Preset {
    name: &'static str,
    unit: &'static str,
    sim: fn() -> Sim,
}

fn three_per_apex() -> Vec<(Vertex, bool, f64)> {
    Vertex::ALL.iter().flat_map(|&v| [0.25, 0.5, 0.75].map(|t| (v, false, t))).collect()
}

/// The two sides a cross cut joins: PQ and PQ′.
fn cross_name(v: Vertex) -> &'static str {
    ["AB×AC", "BC×BA", "CA×CB"][v.index()]
}

const PRESETS: [Preset; 6] = [
    Preset {
        name: "Survey: 100 km triangle on R = 6371 km, σ = 1 cm",
        unit: "km",
        sim: || Sim {
            surface: Surface::Sphere,
            radius: 6371.0,
            tri: Triangle { a: 100.0, b: 90.0, c: 110.0 },
            cuts: three_per_apex(),
            noise: 1e-5,
            seed: 1,
        },
    },
    Preset {
        name: "Gauss 1828 Art. 28: A Hohehagen, B Brocken, C Inselberg (a ≈ 107, b ≈ 85, c ≈ 69 km), R = 6371 km, σ = 1 cm",
        unit: "km",
        sim: || Sim {
            surface: Surface::Sphere,
            radius: 6371.0,
            tri: Triangle { a: 107.0, b: 85.0, c: 69.0 },
            cuts: three_per_apex(),
            noise: 1e-5,
            seed: 1,
        },
    },
    Preset {
        name: "Octant of a sphere, R = 1, exact",
        unit: "u",
        sim: || Sim {
            surface: Surface::Sphere,
            radius: 1.0,
            tri: Triangle { a: PI / 2.0, b: PI / 2.0, c: PI / 2.0 },
            cuts: Vertex::ALL.iter().map(|&v| (v, false, 0.5)).collect(),
            noise: 0.0,
            seed: 1,
        },
    },
    Preset {
        name: "Flat field: 1 km triangle, σ = 1 mm",
        unit: "m",
        sim: || Sim {
            surface: Surface::Plane,
            radius: 1.0,
            tri: Triangle { a: 1000.0, b: 900.0, c: 1100.0 },
            cuts: three_per_apex(),
            noise: 0.001,
            seed: 1,
        },
    },
    Preset {
        name: "Saddle: pseudo-radius 10, sides ≈ 6, σ = 0.001",
        unit: "u",
        sim: || Sim {
            surface: Surface::Saddle,
            radius: 10.0,
            tri: Triangle { a: 6.0, b: 5.4, c: 6.6 },
            cuts: three_per_apex(),
            noise: 0.001,
            seed: 1,
        },
    },
    Preset {
        name: "Survey: 1 km triangle on R = 6371 km, σ = 1 mm",
        unit: "km",
        sim: || Sim {
            surface: Surface::Sphere,
            radius: 6371.0,
            tri: Triangle { a: 1.0, b: 0.9, c: 1.1 },
            cuts: three_per_apex(),
            noise: 1e-6,
            seed: 1,
        },
    },
];

pub struct App {
    unit: String,
    sigma: f64,
    tri: Triangle,
    cuts: Vec<CutInput>,
    sim: Sim,
    /// K of the simulator surface while the measurements are unedited sim output.
    sim_k: Option<f64>,
    sim_msg: String,
    tab: Tab,
    /// Display-only curvature factor for the map and sphere views.
    exaggerate: f64,
    yaw: f32,
    pitch: f32,
    zoom: f32,
    report_msg: String,
}

impl App {
    pub fn new() -> Self {
        let mut app = App {
            unit: "m".into(),
            sigma: 0.0,
            tri: Triangle { a: 0.0, b: 0.0, c: 0.0 },
            cuts: Vec::new(),
            sim: (PRESETS[0].sim)(),
            sim_k: None,
            sim_msg: String::new(),
            tab: Tab::Results,
            exaggerate: 1.0,
            yaw: 0.0,
            pitch: 0.0,
            zoom: 1.0,
            report_msg: String::new(),
        };
        app.load();
        app
    }

    /// Last measurement set, kept across restarts.
    fn save_path() -> Option<std::path::PathBuf> {
        Some(std::path::PathBuf::from(std::env::var_os("HOME")?).join(".local/share/gauss_sim/measurements.txt"))
    }

    fn save(&self) {
        let Some(path) = Self::save_path() else { return };
        let t = &self.tri;
        let mut out = format!("unit {}\nsides {} {} {}\nsigma {}\n", self.unit, t.a, t.b, t.c, self.sigma);
        for c in &self.cuts {
            out += &match c.y {
                None => format!("cut {} {} {}\n", c.apex.name(), c.m, c.d),
                Some(y) => format!("cross {} {} {} {}\n", c.apex.name(), c.m, y, c.d),
            };
        }
        // ponytail: a failed save only loses the convenience of reloading; the session goes on.
        let _ = std::fs::create_dir_all(path.parent().unwrap()).and_then(|_| std::fs::write(&path, out));
    }

    fn load(&mut self) {
        let Some(text) = Self::save_path().and_then(|p| std::fs::read_to_string(p).ok()) else { return };
        for line in text.lines() {
            let w: Vec<&str> = line.split_whitespace().collect();
            let f = |i: usize| w.get(i).and_then(|x| x.parse::<f64>().ok()).unwrap_or(0.0);
            match w.first() {
                Some(&"unit") => self.unit = w.get(1).unwrap_or(&"").to_string(),
                Some(&"sides") => self.tri = Triangle { a: f(1), b: f(2), c: f(3) },
                Some(&"sigma") => self.sigma = f(1),
                Some(&kind @ ("cut" | "cross")) => {
                    let apex = Vertex::ALL.into_iter().find(|v| Some(&v.name()) == w.get(1)).unwrap_or(Vertex::A);
                    let c = if kind == "cut" {
                        CutInput { apex, m: f(2), y: None, d: f(3) }
                    } else {
                        CutInput { apex, m: f(2), y: Some(f(3)), d: f(4) }
                    };
                    self.cuts.push(c);
                }
                _ => {}
            }
        }
    }

    fn load_preset(&mut self, i: usize) {
        self.sim = (PRESETS[i].sim)();
        self.unit = PRESETS[i].unit.into();
        self.generate();
    }

    fn generate(&mut self) {
        match self.sim.measure() {
            Ok((tri, cuts)) => {
                self.tri = tri;
                self.cuts = cuts;
                self.sigma = self.sim.noise;
                self.sim_k = Some(self.sim.k_true());
                self.sim_msg.clear();
                self.reset_view();
                self.save();
            }
            Err(e) => self.sim_msg = e,
        }
    }

    fn reset_view(&mut self) {
        (self.yaw, self.pitch, self.zoom) = (0.0, 0.0, 1.0);
        self.exaggerate = 1.0;
    }
}

/// Measured length: click to type. No drag, so a stray mouse move cannot change a measurement.
fn len_entry(ui: &mut egui::Ui, v: &mut f64) -> egui::Response {
    ui.add(egui::DragValue::new(v).speed(0.0).max_decimals(12).range(0.0..=f64::MAX))
}

fn len_drag(ui: &mut egui::Ui, v: &mut f64) -> egui::Response {
    let speed = (v.abs() * 1e-6).max(1e-12);
    ui.add(egui::DragValue::new(v).speed(speed).max_decimals(12).range(0.0..=f64::MAX))
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        egui::Panel::top("tabs").show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.heading("Gauss_Sim");
                ui.label("curvature from measured distances");
                ui.separator();
                for (t, name) in TABS {
                    ui.selectable_value(&mut self.tab, t, name);
                }
            });
        });

        egui::Panel::left("inputs").resizable(true).default_size(380.0).show(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| self.inputs_ui(ui));
        });

        let an = analyze(&self.tri, &self.cuts, self.sigma);
        egui::CentralPanel::default().show(ui, |ui| match self.tab {
            Tab::Results => egui::ScrollArea::both().show(ui, |ui| self.results_ui(ui, &an)).inner,
            Tab::Ratio => self.ratio_ui(ui, &an),
            Tab::CutFunctions => self.cut_fn_ui(ui, &an),
            Tab::Map => self.map_ui(ui, &an),
            Tab::Sphere => self.sphere_ui(ui, &an),
            Tab::Notation => egui::ScrollArea::vertical().show(ui, notation_ui).inner,
            Tab::Report => self.report_ui(ui, &an),
        });
    }
}

impl App {
    fn inputs_ui(&mut self, ui: &mut egui::Ui) {
        let mut gen = false;
        egui::CollapsingHeader::new("Simulator (makes test data, overwrites Measurements)").default_open(false).show(ui, |ui| {
            egui::ComboBox::from_label("Preset").selected_text("load…").show_ui(ui, |ui| {
                for (i, p) in PRESETS.iter().enumerate() {
                    if ui.selectable_label(false, p.name).clicked() {
                        self.load_preset(i);
                    }
                }
            });
            let s = &mut self.sim;
            ui.horizontal(|ui| {
                ui.radio_value(&mut s.surface, Surface::Plane, "Plane K = 0");
                ui.radio_value(&mut s.surface, Surface::Sphere, "Sphere K > 0");
                ui.radio_value(&mut s.surface, Surface::Saddle, "Saddle K < 0");
            });
            egui::Grid::new("sim").num_columns(2).show(ui, |ui| {
                if s.surface != Surface::Plane {
                    ui.label("R_sim");
                    len_drag(ui, &mut s.radius);
                    ui.end_row();
                }
                for (name, v) in [("a_sim = |BC|", &mut s.tri.a), ("b_sim = |CA|", &mut s.tri.b), ("c_sim = |AB|", &mut s.tri.c)] {
                    ui.label(name);
                    len_drag(ui, v);
                    ui.end_row();
                }
                ui.label("noise σ");
                len_drag(ui, &mut s.noise);
                ui.end_row();
                ui.label("seed");
                ui.add(egui::DragValue::new(&mut s.seed));
                ui.end_row();
            });
            ui.label("Cuts: apex, cross, position t (foot m = tℓ; cross x = tu, y = tv)");
            let mut remove = None;
            for (i, (apex, cross, t)) in s.cuts.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    apex_combo(ui, ("simapex", i), apex);
                    ui.checkbox(cross, "cross");
                    ui.label(if *cross { cross_name(*apex) } else { apex.base_name() });
                    ui.add(egui::Slider::new(t, 0.01..=0.99).text("t"));
                    if ui.small_button("×").clicked() {
                        remove = Some(i);
                    }
                });
            }
            if let Some(i) = remove {
                s.cuts.remove(i);
            }
            ui.horizontal(|ui| {
                if ui.button("+ cut").clicked() {
                    s.cuts.push((Vertex::A, false, 0.5));
                }
                gen = ui.button("Generate measurements ⟶").clicked();
            });
        });
        if gen {
            self.generate();
        }
        if !self.sim_msg.is_empty() {
            ui.colored_label(Color32::LIGHT_RED, &self.sim_msg);
        }

        ui.separator();
        ui.heading("Measurements");
        ui.label("Geodesic lengths measured on the surface. No angle and no radius is entered.");
        let mut edited = false;
        egui::Grid::new("meas").num_columns(2).show(ui, |ui| {
            ui.label("unit");
            edited |= ui.add(egui::TextEdit::singleline(&mut self.unit).desired_width(60.0)).changed();
            ui.end_row();
            for (name, v) in [("a = |BC|", &mut self.tri.a), ("b = |CA|", &mut self.tri.b), ("c = |AB|", &mut self.tri.c)] {
                ui.label(name);
                edited |= len_entry(ui, v).changed();
                ui.end_row();
            }
            ui.label("σ (per length)");
            edited |= len_entry(ui, &mut self.sigma).changed();
            ui.end_row();
        });
        ui.label("Cuts j from apex P.\nfoot: m = |QD| along the base from Q, d = |PD|.\ncross: x = |PX| on PQ, y = |PY| on PQ′, d = |XY|.");
        let mut remove = None;
        egui::Grid::new("cuts").striped(true).num_columns(7).show(ui, |ui| {
            for h in ["j", "apex", "type", "m or x", "y", "d", ""] {
                ui.label(h);
            }
            ui.end_row();
            for (i, c) in self.cuts.iter_mut().enumerate() {
                ui.label(format!("{}", i + 1));
                edited |= apex_combo(ui, ("apex", i), &mut c.apex);
                let mut cross = c.y.is_some();
                let kinds = [(false, format!("foot on {}", c.apex.base_name())), (true, format!("cross {}", cross_name(c.apex)))];
                egui::ComboBox::from_id_salt(("kind", i)).selected_text(&kinds[cross as usize].1).show_ui(ui, |ui| {
                    for (k, name) in &kinds {
                        ui.selectable_value(&mut cross, *k, name);
                    }
                });
                if cross != c.y.is_some() {
                    c.y = cross.then_some(c.m);
                    edited = true;
                }
                edited |= len_entry(ui, &mut c.m).changed();
                match &mut c.y {
                    Some(y) => edited |= len_entry(ui, y).changed(),
                    None => {
                        ui.label("");
                    }
                }
                edited |= len_entry(ui, &mut c.d).changed();
                if ui.small_button("×").clicked() {
                    remove = Some(i);
                }
                ui.end_row();
            }
        });
        if let Some(i) = remove {
            self.cuts.remove(i);
            edited = true;
        }
        if ui.button("+ cut").clicked() {
            let (l, _, _) = self.tri.frame(Vertex::A);
            self.cuts.push(CutInput { apex: Vertex::A, m: l / 2.0, y: None, d: 0.0 });
            edited = true;
        }
        if ui.button("Clear all measurements").clicked() {
            self.tri = Triangle { a: 0.0, b: 0.0, c: 0.0 };
            self.cuts.clear();
            self.sigma = 0.0;
            edited = true;
        }
        if edited {
            self.sim_k = None;
            self.save();
        }
    }

    fn report_ui(&mut self, ui: &mut egui::Ui, an: &Analysis) {
        let text = crate::report::build(&self.tri, &self.cuts, self.sigma, &self.unit, an, self.sim_k);
        ui.horizontal(|ui| {
            if ui.button("Copy report").clicked() {
                ui.ctx().copy_text(text.clone());
                self.report_msg = "Copied to the clipboard.".into();
            }
            if ui.button("Save report as .txt").clicked() {
                let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs());
                let dir = std::path::PathBuf::from(std::env::var_os("HOME").unwrap_or_default()).join("Documents");
                let path = dir.join(format!("gauss_sim_report_{secs}.txt"));
                self.report_msg = match std::fs::create_dir_all(&dir).and_then(|_| std::fs::write(&path, &text)) {
                    Ok(()) => format!("Saved to {}", path.display()),
                    Err(e) => format!("Save failed: {e}"),
                };
            }
            ui.label(&self.report_msg);
        });
        ui.label("Every input, equation, intermediate value and result behind the other tabs. Paste it into any LLM to check the maths.");
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.add(egui::TextEdit::multiline(&mut text.as_str()).font(egui::TextStyle::Monospace).desired_width(f32::INFINITY));
        });
    }

    fn results_ui(&mut self, ui: &mut egui::Ui, an: &Analysis) {
        let u = &self.unit;
        if !an.triangle_ok {
            if self.tri.a == 0.0 && self.tri.b == 0.0 && self.tri.c == 0.0 {
                ui.label("Enter the sides a, b, c and at least one cut in Measurements. Click a box to type a value.");
            } else {
                ui.colored_label(Color32::LIGHT_RED, "The sides a, b, c fail the triangle inequality.");
            }
        }
        let (text, color) = match an.verdict {
            Some(Verdict::Flat) => ("FLAT within measurement: K* is within 2σ(K*) of 0", Color32::from_rgb(220, 220, 120)),
            Some(Verdict::Positive) => ("CURVED, K* > 0: sphere-like (cuts longer than the plane allows)", apex_color(Vertex::A)),
            Some(Verdict::Negative) => ("CURVED, K* < 0: saddle-like (cuts shorter than the plane allows)", apex_color(Vertex::C)),
            None => ("No valid cut: add a cut to fix K", Color32::GRAY),
        };
        ui.label(egui::RichText::new(text).size(20.0).color(color).strong());
        ui.add_space(6.0);

        egui::Grid::new("curv").num_columns(2).striped(true).show(ui, |ui| {
            let Some(k) = an.k_star else { return };
            let sk = an.sigma_k_star;
            ui.label("K*  (combined curvature)");
            ui.label(match sk {
                Some(s) => format!("{} ± {} {u}⁻²", num(k), num(s)),
                None => format!("{} {u}⁻²", num(k)),
            });
            ui.end_row();
            if k != 0.0 {
                let r = 1.0 / k.abs().sqrt();
                ui.label(if k > 0.0 { "R = 1/√K*  (sphere radius)" } else { "R = 1/√(−K*)  (pseudo-radius)" });
                ui.label(match sk {
                    Some(s) => format!("{} ± {} {u}", num(r), num(s / (2.0 * k.abs().powf(1.5)))),
                    None => format!("{} {u}", num(r)),
                });
                ui.end_row();
            }
            ui.label("ω = √K* = 1/R  (turn per unit distance)");
            ui.label(if k > 0.0 {
                let w = k.sqrt();
                format!("{} rad/{u} = {} °/{u} = {} ″/{u}", num(w), num(w.to_degrees()), num(w.to_degrees() * 3600.0))
            } else {
                "not real for K* ≤ 0".into()
            });
            ui.end_row();
            if let (Some(Verdict::Flat), Some(s)) = (an.verdict, sk) {
                ui.label("2σ(K*) limits");
                ui.label(format!(
                    "any sphere that fits: R ≥ {} {u};  any saddle: R ≥ {} {u}",
                    num(1.0 / (k + 2.0 * s).max(0.0).sqrt()),
                    num(1.0 / (2.0 * s - k).max(0.0).sqrt())
                ));
                ui.end_row();
            }
            ui.label("χ²/dof of Kⱼ about K*");
            ui.label(an.chi2_dof.map_or("needs σ > 0 and 2+ cuts".into(), |c| {
                format!("{}  ({})", num(c), if c < 3.0 { "cuts agree on one K" } else { "cuts disagree: K not constant, or a blunder" })
            }));
            ui.end_row();
            if let Some(kt) = self.sim_k {
                ui.label("simulator K_sim");
                let pull = sk.filter(|s| *s > 0.0).map_or(String::new(), |s| format!(",  (K* − K_sim)/σ(K*) = {:+.2}", (k - kt) / s));
                ui.label(format!("K_sim = {} {u}⁻²{}", num(kt), pull));
                ui.end_row();
            }
        });

        ui.add_space(8.0);
        ui.strong("Triangle on K* (derived from lengths, not measured)");
        egui::Grid::new("tri").num_columns(2).striped(true).show(ui, |ui| {
            if let (Some([al, be, ga]), Some(area)) = (an.angles, an.area) {
                let k = an.k_star.unwrap_or(0.0);
                for (name, x) in [("α at A", al), ("β at B", be), ("γ at C", ga)] {
                    ui.label(name);
                    ui.label(format!("{:.9}°", x.to_degrees()));
                    ui.end_row();
                }
                ui.label("α + β + γ");
                ui.label(format!("{:.9}°", (al + be + ga).to_degrees()));
                ui.end_row();
                let eps = k * area;
                ui.label("ε = K*𝒜  (angle excess, Gauss–Bonnet)");
                ui.label(format!("{} rad = {} ″", num(eps), num(eps.to_degrees() * 3600.0)));
                ui.end_row();
                ui.label("𝒜  (area)");
                ui.label(format!("{} {u}²", num(area)));
                ui.end_row();
            }
        });

        ui.add_space(8.0);
        ui.strong("Per cut");
        egui::Grid::new("percut").striped(true).show(ui, |ui| {
            for h in ["j", "P", "type", "m | x", "n | y", "d", "d₀ (flat)", "ρ − 1", "Δ", "Kⱼ", "σ(Kⱼ)", "Rⱼ", "d − d(K*)"] {
                ui.strong(h);
            }
            ui.end_row();
            for (j, c) in an.cuts.iter().enumerate() {
                let cut = &c.cut;
                ui.label(format!("{}", j + 1));
                ui.colored_label(apex_color(c.apex), c.apex.name());
                ui.label(if cut.y.is_some() { format!("cross {}", cross_name(c.apex)) } else { format!("foot {}", c.apex.base_name()) });
                ui.label(num(cut.m));
                ui.label(num(cut.y.unwrap_or(cut.n())));
                ui.label(num(cut.d));
                ui.label(num(cut.d_flat()));
                ui.label(num(cut.rho() - 1.0));
                ui.label(num(cut.defect()));
                match c.k {
                    Ok(k) => {
                        ui.label(num(k));
                        ui.label(num(c.sigma_k));
                        ui.label(if k == 0.0 { "∞".into() } else { num(1.0 / k.abs().sqrt()) });
                        ui.label(an.k_star.map_or("n/a".into(), |kh| num(cut.d - cut.length_on(kh))));
                    }
                    Err(e) => {
                        ui.colored_label(Color32::LIGHT_RED, e.to_string());
                    }
                }
                ui.end_row();
            }
        });
    }

    fn ratio_ui(&mut self, ui: &mut egui::Ui, an: &Analysis) {
        let max_dev = an.cuts.iter().map(|c| (c.cut.rho() - 1.0).abs()).filter(|x| x.is_finite()).fold(0.0, f64::max);
        let (scale, label) = if max_dev < 1e-3 { (1e6, "ρ − 1  (ppm)") } else { (100.0, "ρ − 1  (%)") };
        ui.label("ρ = d / d₀: measured cut over the flat-plane (Stewart) value d₀. Flat: ρ = 1 everywhere. Bars: ±2σ.");
        Plot::new("ratio")
            .legend(Legend::default())
            .x_axis_label("t = m / ℓ (foot) or x / u (cross)")
            .y_axis_label(label)
            .show(ui, |p| {
                p.hline(HLine::new("flat: ρ = 1", 0.0).color(Color32::GRAY).style(LineStyle::dashed_loose()));
                for (v, cross) in Vertex::ALL.into_iter().flat_map(|v| [(v, false), (v, true)]) {
                    let mine: Vec<_> = an.cuts.iter().filter(|c| c.apex == v && c.cut.y.is_some() == cross).collect();
                    if mine.is_empty() {
                        continue;
                    }
                    let name = if cross { format!("cross {}", cross_name(v)) } else { format!("foot from {}", v.name()) };
                    let (l, u_, v_) = self.tri.frame(v);
                    // Cross-cut curve: x = tu, y = tv.
                    let curve = |k: f64| -> PlotPoints<'static> {
                        (1..100)
                            .map(|i| {
                                let t = i as f64 / 100.0;
                                let (m, y) = if cross { (t * u_, Some(t * v_)) } else { (t * l, None) };
                                let c = Cut { base: l, near: u_, far: v_, m, y, d: 0.0 };
                                [t, (c.length_on(k) / c.d_flat() - 1.0) * scale]
                            })
                            .collect()
                    };
                    let style = if cross { LineStyle::dashed_dense() } else { LineStyle::Solid };
                    if let Some(k) = an.k_star {
                        p.line(Line::new(format!("model K*, {name}"), curve(k)).color(apex_color(v)).style(style));
                    }
                    if let Some(kt) = self.sim_k {
                        p.line(Line::new(format!("K_sim, {name}"), curve(kt)).color(apex_color(v)).style(LineStyle::dotted_dense()));
                    }
                    let pts: Vec<[f64; 2]> = mine.iter().map(|c| [c.cut.t(), (c.cut.rho() - 1.0) * scale]).collect();
                    for c in &mine {
                        let (x, y, e) = (c.cut.t(), (c.cut.rho() - 1.0) * scale, 2.0 * c.sigma_rho * scale);
                        if e > 0.0 {
                            p.line(Line::new(name.clone(), vec![[x, y - e], [x, y + e]]).color(apex_color(v)));
                        }
                    }
                    let shape = if cross { egui_plot::MarkerShape::Diamond } else { egui_plot::MarkerShape::Circle };
                    p.points(Points::new(name, pts).radius(5.0).shape(shape).color(apex_color(v)));
                }
            });
    }

    fn cut_fn_ui(&mut self, ui: &mut egui::Ui, an: &Analysis) {
        ui.label("Gₖ(K) for each cut (foot cuts scaled by ℓ³). Each curve crosses zero once, at that cut's Kⱼ. Gₖ > 0 left of the root, < 0 right of it.");
        let ks: Vec<f64> = an.cuts.iter().filter_map(|c| c.k.ok()).collect();
        let Some(kh) = an.k_star else {
            ui.label("No valid cut.");
            return;
        };
        let spread = ks.iter().map(|k| (k - kh).abs()).fold(0.0, f64::max);
        let half = [3.0 * an.sigma_k_star.unwrap_or(0.0), 1.5 * kh.abs(), 2.0 * spread, 1e-9 / self.tri.scale().powi(2)]
            .into_iter()
            .fold(0.0, f64::max);
        Plot::new("gk").legend(Legend::default()).x_axis_label(format!("K  ({}⁻²)", self.unit)).show(ui, |p| {
            p.hline(HLine::new("Gₖ = 0", 0.0).color(Color32::GRAY));
            for (j, c) in an.cuts.iter().enumerate() {
                if c.k.is_err() {
                    continue;
                }
                let kmax = c.cut.k_max();
                let l3 = if c.cut.y.is_some() { 1.0 } else { c.cut.base.powi(3) };
                let pts: PlotPoints = (0..=200)
                    .map(|i| kh - half + 2.0 * half * i as f64 / 200.0)
                    .filter(|&k| k < kmax)
                    .map(|k| [k, c.cut.g(k) / l3])
                    .collect();
                p.line(Line::new(format!("j = {} (apex {})", j + 1, c.apex.name()), pts).color(apex_color(c.apex)));
            }
        });
    }

    /// Display curvature: K times the exaggeration, capped so the triangle still fits a sphere.
    fn k_draw(&self, an: &Analysis) -> f64 {
        let k = an.k_star.unwrap_or(0.0) * self.exaggerate;
        let p = self.tri.a + self.tri.b + self.tri.c;
        // Saddle side: past |K|·p² ≈ 1600 the hyperboloid coordinates lose the triangle.
        k.clamp(-(40.0 / p).powi(2), (1.9 * PI / p).powi(2))
    }

    fn exaggerate_ui(&mut self, ui: &mut egui::Ui, an: &Analysis) {
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut self.exaggerate, 1.0..=1e6).logarithmic(true).text("display × K*"));
            if self.exaggerate != 1.0 {
                ui.colored_label(Color32::YELLOW, format!("display only: drawn with K = {} (not the measured K*)", num(self.k_draw(an))));
            }
        });
    }

    fn map_ui(&mut self, ui: &mut egui::Ui, an: &Analysis) {
        if an.angles.is_none() {
            ui.label("No valid triangle.");
            return;
        }
        self.exaggerate_ui(ui, an);
        ui.label("Azimuthal-equidistant map about A: every distance from A and every angle at A is true. On a plane every geodesic is straight; the gap between side BC and the dashed chord B–C shows the curvature.");
        let model = Model { k: self.k_draw(an) };
        let geo = Geo::new(&model, &self.tri, &self.cuts);
        let flat = |p: P3| {
            let (r, az) = model.polar(p);
            [r * az.cos(), r * az.sin()]
        };
        let paths = geo.paths(&model);
        let mapped: Vec<(Color32, f32, Vec<[f64; 2]>)> = paths.iter().map(|(c, w, pts)| (*c, *w, pts.iter().map(|&p| flat(p)).collect())).collect();

        let (resp, painter) = ui.allocate_painter(ui.available_size(), egui::Sense::hover());
        let rect = resp.rect.shrink(30.0);
        let all = mapped.iter().flat_map(|m| m.2.iter());
        let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
        for p in all {
            for i in 0..2 {
                lo[i] = lo[i].min(p[i]);
                hi[i] = hi[i].max(p[i]);
            }
        }
        let s = (rect.width() as f64 / (hi[0] - lo[0])).min(rect.height() as f64 / (hi[1] - lo[1]));
        let to = |p: [f64; 2]| {
            Pos2::new(
                rect.center().x + ((p[0] - (lo[0] + hi[0]) / 2.0) * s) as f32,
                rect.center().y - ((p[1] - (lo[1] + hi[1]) / 2.0) * s) as f32,
            )
        };
        // Flat chord B–C for reference.
        let (b2, c2) = (to(flat(geo.vb)), to(flat(geo.vc)));
        painter.extend(egui::Shape::dashed_line(&[b2, c2], Stroke::new(1.0, Color32::GRAY), 6.0, 4.0));
        for (col, w, pts) in &mapped {
            painter.add(egui::Shape::line(pts.iter().map(|&p| to(p)).collect(), Stroke::new(*w, *col)));
        }
        for (name, p) in [("A", geo.va), ("B", geo.vb), ("C", geo.vc)] {
            let q = to(flat(p));
            painter.circle_filled(q, 4.0, Color32::WHITE);
            painter.text(q + Vec2::new(8.0, -8.0), egui::Align2::LEFT_BOTTOM, name, egui::FontId::proportional(18.0), Color32::WHITE);
        }
        for (j, (apex, _, foot)) in geo.feet.iter().enumerate() {
            let q = to(flat(*foot));
            painter.circle_filled(q, 3.0, apex_color(*apex));
            painter.text(q + Vec2::new(5.0, 5.0), egui::Align2::LEFT_TOP, format!("D{}", subscript(j + 1)), egui::FontId::proportional(13.0), apex_color(*apex));
        }
    }

    fn sphere_ui(&mut self, ui: &mut egui::Ui, an: &Analysis) {
        if an.angles.is_none() {
            ui.label("No valid triangle.");
            return;
        }
        self.exaggerate_ui(ui, an);
        let k = self.k_draw(an);
        if k <= 0.0 {
            ui.label(if k == 0.0 {
                "K* = 0: the model surface is the plane (see Intrinsic map)."
            } else {
                "K* < 0: a complete surface of constant negative curvature has no smooth isometric embedding in 3-D space (Hilbert, 1901). See Intrinsic map."
            });
            return;
        }
        ui.horizontal(|ui| {
            ui.label("Drag to rotate, scroll to zoom.");
            if ui.button("Reset view").clicked() {
                (self.yaw, self.pitch, self.zoom) = (0.0, 0.0, 1.0);
            }
            ui.label(format!("R drawn = {} {}", num(1.0 / k.sqrt()), self.unit));
        });
        let model = Model { k };
        let r = 1.0 / k.sqrt();
        let geo = Geo::new(&model, &self.tri, &self.cuts);

        let (resp, painter) = ui.allocate_painter(ui.available_size(), egui::Sense::drag());
        let d = resp.drag_delta();
        self.yaw += d.x * 0.01;
        self.pitch = (self.pitch + d.y * 0.01).clamp(-1.5, 1.5);
        if resp.hovered() {
            let sc = ui.input(|i| i.smooth_scroll_delta.y);
            self.zoom = (self.zoom * (sc * 0.003).exp()).clamp(0.05, 1e7);
        }

        // View basis: centroid of the triangle faces the viewer.
        let unit = |p: P3| p.map(|x| x / r);
        let g = norm(add(add(unit(geo.va), unit(geo.vb)), unit(geo.vc)));
        let e1 = norm(cross(if g[2].abs() < 0.9 { [0.0, 0.0, 1.0] } else { [1.0, 0.0, 0.0] }, g));
        let e2 = cross(g, e1);
        let (sy, cy) = (self.yaw as f64).sin_cos();
        let (sp, cp) = (self.pitch as f64).sin_cos();
        let view = |p: P3| {
            let v = [dot(p, e1), dot(p, e2), dot(p, g)];
            let v = [cy * v[0] + sy * v[2], v[1], -sy * v[0] + cy * v[2]];
            [v[0], cp * v[1] - sp * v[2], sp * v[1] + cp * v[2]]
        };
        // Fit the triangle: its spread on the unit sphere sets the default zoom.
        let spread = [geo.va, geo.vb, geo.vc].iter().map(|&p| norm2(sub(unit(p), g))).fold(1e-12, f64::max).sqrt();
        let rect = resp.rect;
        let base = (rect.width().min(rect.height()) * 0.45) as f64 * (0.8 / spread).max(1.0) * self.zoom as f64;
        let to = |v: [f64; 3]| Pos2::new(rect.center().x + (v[0] * base) as f32, rect.center().y - (v[1] * base) as f32);

        painter.circle_stroke(rect.center(), base as f32, Stroke::new(1.0, Color32::from_gray(90)));
        let grid = Color32::from_gray(60);
        let draw = |pts: &[P3], col: Color32, w: f32| {
            for s in pts.windows(2) {
                let (a, b) = (view(unit(s[0])), view(unit(s[1])));
                let front = a[2] >= 0.0 && b[2] >= 0.0;
                let c = if front { col } else { col.gamma_multiply(0.25) };
                painter.line_segment([to(a), to(b)], Stroke::new(if front { w } else { 1.0 }, c));
            }
        };
        for i in 1..12 {
            let lat = model_lat_circle(&model, i as f64 * PI / 12.0 * r);
            draw(&lat, grid, 1.0);
        }
        for i in 0..24 {
            let az = i as f64 * PI / 12.0;
            let mer: Vec<P3> = (0..=48).map(|j| model.polar_point(j as f64 / 48.0 * PI * r, az)).collect();
            draw(&mer, grid, 1.0);
        }
        for (col, w, pts) in geo.paths(&model) {
            draw(&pts, col, w);
        }
        for (name, p) in [("A", geo.va), ("B", geo.vb), ("C", geo.vc)] {
            let v = view(unit(p));
            if v[2] >= 0.0 {
                painter.circle_filled(to(v), 4.0, Color32::WHITE);
                painter.text(to(v) + Vec2::new(8.0, -8.0), egui::Align2::LEFT_BOTTOM, name, egui::FontId::proportional(18.0), Color32::WHITE);
            }
        }
    }
}

fn model_lat_circle(model: &Model, dist: f64) -> Vec<P3> {
    (0..=72).map(|j| model.polar_point(dist, j as f64 / 72.0 * 2.0 * PI)).collect()
}

fn subscript(n: usize) -> String {
    n.to_string().chars().map(|c| char::from_u32(0x2080 + c.to_digit(10).unwrap()).unwrap()).collect()
}

fn add(a: P3, b: P3) -> P3 {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: P3, b: P3) -> P3 {
    std::array::from_fn(|i| a[i] - b[i])
}
fn dot(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm2(a: P3) -> f64 {
    dot(a, a)
}
fn norm(a: P3) -> P3 {
    let l = norm2(a).sqrt();
    a.map(|x| x / l)
}
fn cross(a: P3, b: P3) -> P3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn apex_combo(ui: &mut egui::Ui, id: impl std::hash::Hash + std::fmt::Debug, v: &mut Vertex) -> bool {
    let before = *v;
    egui::ComboBox::from_id_salt(id).width(40.0).selected_text(v.name()).show_ui(ui, |ui| {
        for x in Vertex::ALL {
            ui.selectable_value(v, x, x.name());
        }
    });
    before != *v
}

/// The measured triangle placed in the model surface: A at the pole, B on azimuth 0, C on azimuth α.
struct Geo {
    va: P3,
    vb: P3,
    vc: P3,
    /// (apex, start, end) per cut: P → Dⱼ, or X → Y for a cross cut.
    feet: Vec<(Vertex, P3, P3)>,
}

impl Geo {
    fn new(model: &Model, tri: &Triangle, cuts: &[CutInput]) -> Self {
        let alpha = crate::geom::angle(model.k, tri.b, tri.c, tri.a);
        let va = model.polar_point(0.0, 0.0);
        let vb = model.polar_point(tri.c, 0.0);
        let vc = model.polar_point(tri.b, alpha);
        let v = [va, vb, vc];
        let feet = cuts
            .iter()
            .map(|c| {
                let i = c.apex.index();
                let (q, q2) = (v[(i + 1) % 3], v[(i + 2) % 3]);
                let (l, u, w) = tri.frame(c.apex);
                match c.y {
                    None => (c.apex, v[i], model.lerp(q, q2, (c.m / l).clamp(0.0, 1.0))),
                    Some(y) => (c.apex, model.lerp(v[i], q, (c.m / u).clamp(0.0, 1.0)), model.lerp(v[i], q2, (y / w).clamp(0.0, 1.0))),
                }
            })
            .collect();
        Geo { va, vb, vc, feet }
    }

    fn paths(&self, model: &Model) -> Vec<(Color32, f32, Vec<P3>)> {
        let seg = |p: P3, q: P3| (0..=64).map(|i| model.lerp(p, q, i as f64 / 64.0)).collect::<Vec<_>>();
        let mut out = vec![
            (Color32::WHITE, 2.5, seg(self.va, self.vb)),
            (Color32::WHITE, 2.5, seg(self.vb, self.vc)),
            (Color32::WHITE, 2.5, seg(self.vc, self.va)),
        ];
        for (apex, start, end) in &self.feet {
            out.push((apex_color(*apex), 1.5, seg(*start, *end)));
        }
        out
    }
}

const NOTATION: &[(&str, &str)] = &[
    ("A, B, C", "vertices of the triangle (marks on the surface)"),
    ("a, b, c", "measured sides: a = |BC|, b = |CA|, c = |AB|"),
    ("Pⱼ", "apex of cut j (A, B or C)"),
    ("(P, Q, Q′)", "cyclic order: (A, B, C), (B, C, A), (C, A, B). The base of cut j is QQ′, the side opposite P"),
    ("foot cut", "from the apex P to a mark Dⱼ on the base QQ′"),
    ("cross cut", "from a mark X on PQ to a mark Y on PQ′ (both sides that meet at P)"),
    ("Dⱼ", "foot of cut j, a mark on the base"),
    ("ℓⱼ = |QQ′|", "base length (a, b or c)"),
    ("uⱼ = |PQ|, vⱼ = |PQ′|", "the two other sides, from the apex"),
    ("mⱼ = |QDⱼ|", "measured along the base, from Q"),
    ("nⱼ = ℓⱼ − mⱼ", "rest of the base, |DⱼQ′|"),
    ("dⱼ = |PⱼDⱼ|", "measured cut length"),
    ("xⱼ = |PX|, yⱼ = |PY|", "cross cut: marks measured from P along PQ and PQ′"),
    ("dⱼ = |XY|", "cross cut length"),
    ("tⱼ", "position: m/ℓ (foot) or x/u (cross)"),
    ("σ", "standard deviation of one length measurement (same for every length)"),
    ("K", "Gaussian curvature, unit⁻². K > 0 sphere-like, K = 0 flat, K < 0 saddle-like"),
    ("Kⱼ", "curvature from cut j alone"),
    ("K*", "combined curvature: inverse-variance mean of the Kⱼ"),
    ("σ(Kⱼ), σ(K*)", "standard deviations from σ; σ(K*) uses every length once, so shared sides count once"),
    ("R = 1/√|K*|", "radius of curvature (sphere radius, or saddle pseudo-radius)"),
    ("Rⱼ = 1/√|Kⱼ|", "radius from cut j alone"),
    ("ω = √K* = 1/R", "turn per unit distance, rad/unit (K* > 0)"),
    ("K_sim, R_sim, a_sim, b_sim, c_sim", "simulator surface and true sides (before noise)"),
    ("subscript ₖ", "evaluated at curvature K: snₖ, csₖ, Gₖ"),
    ("snₖ(x), csₖ(x)", "sin(√K x)/√K, cos(√K x);  x, 1;  sinh(√−K x)/√−K, cosh(√−K x)"),
    ("d₀", "cut length that a plane (K = 0) predicts from the other lengths"),
    ("Δ = d² − d₀²", "defect, unit². Zero on a plane"),
    ("ρ = d / d₀", "ratio test. ρ = 1 flat, ρ > 1 for K > 0, ρ < 1 for K < 0"),
    ("d(K)", "cut length that (2) predicts on curvature K; d₀ = d(0)"),
    ("Gₖ", "cut function; its root is Kⱼ"),
    ("α, β, γ", "angles at A, B, C. Outputs only, computed from lengths and K*"),
    ("ε = α + β + γ − π", "angle excess, ε = K𝒜"),
    ("𝒜", "triangle area, unit²"),
];

const EQUATIONS: &[(&str, &str)] = &[
    ("(1) law of cosines", "csₖ(a) = csₖ(b)·csₖ(c) + K·snₖ(b)·snₖ(c)·cos α"),
    ("(2) cut relation", "csₖ(v)·snₖ(m) + csₖ(u)·snₖ(n) = csₖ(d)·snₖ(ℓ)"),
    ("(2x) cross relation", "[csₖ(ℓ) − csₖ(u)csₖ(v)] / [K snₖ(u)snₖ(v)] = [csₖ(d) − csₖ(x)csₖ(y)] / [K snₖ(x)snₖ(y)]  (both = cos of the angle at P)"),
    ("(3) flat limit (Stewart)", "v²m + u²n = ℓ(d² + mn)   ⇒   d₀² = (v²m + u²n)/ℓ − mn"),
    ("(3x) flat cross cut", "d₀² = x² + y² − xy(u² + v² − ℓ²)/(uv)"),
    ("(4) defect", "Δ = d² − d₀²"),
    ("(5) cut function", "Gₖ = [csₖ(v)snₖ(m) + csₖ(u)snₖ(n) − csₖ(d)snₖ(ℓ)] / K,   G₀ = ℓΔ/2"),
    ("(5x) cross function", "Gₖ = left side of (2x) − right side of (2x)"),
    ("(6) solve", "Gₖ(Kⱼ) = 0,   sign Kⱼ = sign Δⱼ = sign(ρⱼ − 1)"),
    ("(7) Gauss–Bonnet", "ε = α + β + γ − π = K𝒜"),
];

fn notation_ui(ui: &mut egui::Ui) {
    ui.heading("Method: cut the triangle, cancel the angle");
    for line in [
        "1. Measure the three sides a, b, c on the surface. Three lengths alone do not fix K: a triangle with these sides exists on every surface of constant K (up to a size limit for K > 0).",
        "2. Put a mark D on one side. Measure m along that side from its start Q, and measure the cut d from the opposite vertex P to D.",
        "3. The cut makes a sub-triangle PQD. It shares the angle at Q with the full triangle PQQ′. Write the law of cosines (1) for both triangles and eliminate that angle. The result (2) holds only lengths and K.",
        "3x. Cross cut: join a mark X on PQ to a mark Y on PQ′. The small triangle PXY shares the angle at P with PQQ′. Eliminating it gives (2x). Marks at the midpoints of both sides give d₀ = ℓ/2 on a plane.",
        "4. Solve (2) or (2x) for K. Each cut gives one Kⱼ. More cuts (other positions, other apexes) give more values, combined into K*. If they agree (χ²/dof ≈ 1), one curvature fits the whole triangle.",
        "5. On a plane, (2) becomes Stewart's theorem (3). The ratio ρ = d/d₀ is the direct test: ρ = 1 is flat, ρ > 1 is sphere-like, ρ < 1 is saddle-like.",
        "6. From K*: R = 1/√|K*|, and ω = 1/R is the turn in radians per unit distance. The angles α, β, γ and the excess ε are outputs.",
        "Gauss (1828): K depends only on lengths measured on the surface (Theorema Egregium). This program uses only such lengths.",
    ] {
        ui.label(line);
        ui.add_space(2.0);
    }
    ui.add_space(8.0);
    ui.heading("Notation");
    egui::Grid::new("notation").striped(true).num_columns(2).show(ui, |ui| {
        for (s, m) in NOTATION {
            ui.monospace(*s);
            ui.label(*m);
            ui.end_row();
        }
    });
    ui.add_space(8.0);
    ui.heading("Equations");
    egui::Grid::new("eqs").striped(true).num_columns(2).show(ui, |ui| {
        for (n, e) in EQUATIONS {
            ui.label(*n);
            ui.monospace(*e);
            ui.end_row();
        }
    });
}
