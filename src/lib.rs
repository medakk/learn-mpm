extern crate nalgebra_glm as glm;

#[derive(Debug)]
pub struct Particle {
    pub x: glm::Vec2,
    pub v: glm::Vec2,
    pub mass: f32,
}

impl Particle {
    pub fn draw_to_buffer(&self, buffer: &mut Vec<u32>, w: usize, h: usize) {
        let r: i32 = 4;
        for i in -r..r+1 {
            for j in -r..r+1 {
                let x = self.x.x as i32 + i;
                let y = self.x.y as i32 + j;
                if x >= w as i32 || x < 0 || y >= h as i32 || y < 0 {
                    continue;
                }

                let idx = y as usize *w + x as usize;
                buffer[idx] = 0xffff0000;
            }
        }
    }
}

pub struct MPM {
    pub particles: Vec<Particle>,
}

impl MPM {
    pub fn new() -> MPM {
        let mut particles = Vec::new();
        for i in 0..10 {
            let p = Particle{
                x: glm::vec2(i as f32 * 640.0 / 10.0, 480.0 / 2.0),
                v: glm::vec2(0.0, 0.0),
                mass: 1.0,
            };
            println!("{:?}", p);
            particles.push(p);
        }

        MPM {
            particles,
        }
    }
}