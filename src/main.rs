extern crate minifb;

use minifb::{Key, Window, WindowOptions};
use learn_mpm::MPM;

const WIDTH: usize = 640;
const HEIGHT: usize = 640;

fn main() {
    let mut buffer: Vec<u32> = vec![0; WIDTH * HEIGHT];

    let mut window = Window::new(
        "MPM",
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
        let n_iters = 4;
        for _ in 0..n_iters {
            mpm.update(0.1);
        }

        // Clear the buffer
        for i in buffer.iter_mut() {
            *i = 0;
        }

        // draw
        for p in mpm.particles.iter() {
            p.draw_to_buffer(&mut buffer, WIDTH, HEIGHT);
        }

        // check input
        window.get_keys().map(|keys| {
            for t in keys {
                match t {
                    Key::R => { mpm.reset(); },
                    _ => (),
                }
            }
        });

        window
            .update_with_buffer(&buffer, WIDTH, HEIGHT)
            .unwrap();
    }
}