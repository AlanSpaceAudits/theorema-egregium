mod app;
mod geom;
mod report;
mod survey;

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default().with_inner_size([1400.0, 900.0]).with_title("Gauss_Sim"),
        ..Default::default()
    };
    eframe::run_native("Gauss_Sim", options, Box::new(|_cc| Ok(Box::new(app::App::new()))))
}
