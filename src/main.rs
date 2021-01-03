extern crate minifb;
extern crate nalgebra_glm as glm;

use minifb::{Key, Window, WindowOptions};
use learn_mpm::MPM;

const WIDTH: usize = 640;
const HEIGHT: usize = 360;

fn main() {
    let mut buffer: Vec<u32> = vec![0; WIDTH * HEIGHT];

    let mut window = Window::new(
        "Test - ESC to exit",
        WIDTH,
        HEIGHT,
        WindowOptions::default(),
    )
    .unwrap_or_else(|e| {
        panic!("{}", e);
    });

    let mut mpm = MPM::new();

    // Limit to max ~60 fps update rate
    window.limit_update_rate(Some(std::time::Duration::from_micros(16600)));

    while window.is_open() && !window.is_key_down(Key::Escape) {
        // Clear the buffer
        for i in buffer.iter_mut() {
            *i = 0;
        }

        for p in mpm.particles.iter() {
            p.draw_to_buffer(&mut buffer, WIDTH, HEIGHT);
        }

        // We unwrap here as we want this code to exit if it fails. Real applications may want to handle this in a different way
        window
            .update_with_buffer(&buffer, WIDTH, HEIGHT)
            .unwrap();
    }
}