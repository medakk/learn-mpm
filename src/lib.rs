#![allow(non_snake_case)]

extern crate rand;
extern crate nalgebra_glm as glm;

const CELL_DIM: usize = 64;
const GRAVITY: f32 = 0.3;
const ELASTIC_LAMBDA: f32 = 10.0;
const ELASTIC_MU: f32 = 20.0;

fn to_uint_vec(v: &glm::Vec2) -> glm::UVec2 {
    return glm::vec2(v.x as u32, v.y as u32);
}
fn to_f32_vec(v: &glm::UVec2) -> glm::Vec2 {
    return glm::vec2(v.x as f32, v.y as f32);
}

#[derive(Debug, Clone)]
pub struct Particle {
    pub x: glm::Vec2,
    pub v: glm::Vec2,
    pub C: glm::Mat2,
    pub F: glm::Mat2,
    pub volume_0: f32,
    pub mass: f32,
}

impl Particle {
    pub fn draw_to_buffer(&self, buffer: &mut Vec<u32>, w: usize, h: usize) {
        let x = self.x.x * w as f32 / CELL_DIM as f32;
        let y = self.x.y * h as f32 / CELL_DIM as f32;
        let r: i32 = 1;
        for i in -r..r+1 {
            for j in -r..r+1 {
                let x = x as i32 + i;
                let y = y as i32 + j;
                if x >= w as i32 || x < 0 || y >= h as i32 || y < 0 {
                    continue;
                }

                let idx = y as usize *w + x as usize;
                buffer[idx] = 0xffff0000;
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct Cell {
    pub v: glm::Vec2,
    pub mass: f32,
}

impl Cell {
    fn new() -> Cell {
        Cell {
            v: glm::vec2(0.0, 0.0),
            mass: 0.0,
        }
    }

    fn reset(&mut self) {
        self.v = glm::vec2(0.0, 0.0);
        self.mass = 0.0;
    }
}

pub struct MPM {
    pub particles: Vec<Particle>,
    pub cells: Vec<Cell>,
}

impl MPM {
    pub fn new() -> MPM {
        let mut mpm = MPM {
            particles: Vec::new(),
            cells: Vec::new(),
        };
        mpm.reset();

        mpm
    }

    pub fn reset(&mut self) {
        let cells = vec![Cell::new(); CELL_DIM*CELL_DIM];
        let mut particles = Vec::new();
        let n = 100;
        for i in 0..n {
            for j in 0..n {
                let x = CELL_DIM as f32 * glm::vec2(
                    (i as f32 / n as f32)*0.6+0.1,
                    (j as f32 / n as f32)*0.6+0.1);
                let p = Particle {
                    x: x,
                    v: 2.0*glm::vec2(rand::random::<f32>()-0.5, rand::random::<f32>()-0.5),
                    C: glm::zero(),
                    F: glm::mat2(1.0, 0.0, 0.0, 1.0),
                    volume_0: 0.0,
                    mass: 1.0,
                };
                particles.push(p);
            }
        }

        self.cells = cells;
        self.particles = particles;

        self.compute_volume_0();
    }

    pub fn update(&mut self, dt: f32) {
        self.reset_cells();

        self.p2g(dt);
        self.cells_update(dt);
        self.g2p(dt);
    }

    fn compute_volume_0(&mut self) {
        self.p2g(0.0);

        for p in self.particles.iter_mut() {
            let cell_idx: glm::UVec2 = to_uint_vec(&p.x);
            let cell_diff: glm::Vec2 = (&p.x - to_f32_vec(&cell_idx)) - glm::vec2(0.5, 0.5);

            let weights = [
                0.5 * glm::pow(&(glm::vec2(0.5, 0.5) - &cell_diff), &glm::vec2(2.0, 2.0)),
                glm::vec2(0.75, 0.75) - glm::pow(&cell_diff, &glm::vec2(2.0, 2.0)),
                0.5 * glm::pow(&(glm::vec2(0.5, 0.5) + &cell_diff), &glm::vec2(2.0, 2.0)),
            ];

            let mut density = 0.0;
            for gx in 0..3 {
                for gy in 0..3 {
                    let weight = weights[gx].x * weights[gy].y;
                    let cell_x = glm::vec2(cell_idx.x + gx as u32 - 1,
                                           cell_idx.y + gy as u32 - 1);
                    let cell_idx = cell_x.x as usize * CELL_DIM + cell_x.y as usize;
                    if cell_idx >= self.cells.len() {
                        continue;
                    }
                    density += self.cells[cell_idx].mass * weight;
                }
            }

            p.volume_0 = p.mass / density;
        }
    }

    fn reset_cells(&mut self) {
        for cell in self.cells.iter_mut() {
            cell.reset();
        }
    }

    fn p2g(&mut self, dt: f32) {
        for p in self.particles.iter_mut() {
            let F = &p.F;
            let J = glm::determinant(F);
            let volume = p.volume_0 * J;

            // convenience
            let F_T = glm::transpose(F);
            let F_inv_T = glm::inverse(&F_T);
            let F_minus_F_inv_T = F - &F_inv_T;

            // mpm equations
            let P_term_0 = ELASTIC_MU * &F_minus_F_inv_T;
            let P_term_1 = ELASTIC_LAMBDA * J.ln() * &F_inv_T;
            let P = P_term_0 + P_term_1;

            let stress = (1.0 / J) * &P * &F_T;

            let eq_16_term_0 = -volume * 4.0 * stress * dt;

            let cell_idx: glm::UVec2 = to_uint_vec(&p.x);
            let cell_diff: glm::Vec2 = (&p.x - to_f32_vec(&cell_idx)) - glm::vec2(0.5, 0.5);

            let weights = [
                0.5 * glm::pow(&(glm::vec2(0.5, 0.5) - &cell_diff), &glm::vec2(2.0, 2.0)),
                glm::vec2(0.75, 0.75) - glm::pow(&cell_diff, &glm::vec2(2.0, 2.0)),
                0.5 * glm::pow(&(glm::vec2(0.5, 0.5) + &cell_diff), &glm::vec2(2.0, 2.0)),
            ];

            for gx in 0..3 {
                for gy in 0..3 {
                    let weight = weights[gx].x * weights[gy].y;
                    let cell_x = glm::vec2(cell_idx.x + gx as u32 -1,
                                           cell_idx.y + gy as u32 - 1);
                    let cell_dist = (to_f32_vec(&cell_x) - &p.x) + glm::vec2(0.5, 0.5);
                    let Q = &p.C * &cell_dist;

                    let mass_contrib = weight * p.mass;

                    let cell_idx = cell_x.x as usize * CELL_DIM + cell_x.y as usize;
                    if cell_idx >= self.cells.len() {
                        continue;
                    }

                    let cell = &mut self.cells[cell_idx];
                    cell.mass += mass_contrib;
                    cell.v += mass_contrib * (&p.v + Q);

                    let momentum = weight * &eq_16_term_0 * &cell_dist;
                    cell.v += momentum;
                }
            }
        }
    }

    fn cells_update(&mut self, dt: f32) {
        for (i, cell) in self.cells.iter_mut().enumerate() {
            if cell.mass <= 0.0 {
                continue;
            }

            cell.v /= cell.mass;
            cell.v += dt * glm::vec2(0.0, GRAVITY);

            let x = i / CELL_DIM;
            let y = i % CELL_DIM;
            if x < 2 || x > CELL_DIM - 3  {
                cell.v.x = 0.0;
            }
            if y < 2 || y > CELL_DIM - 3  {
                cell.v.y = 0.0;
            }
        }
    }

    fn g2p(&mut self, dt: f32) {
        for p in self.particles.iter_mut() {
            p.v = glm::zero();

            let cell_idx: glm::UVec2 = to_uint_vec(&p.x);
            let cell_diff: glm::Vec2 = (&p.x - to_f32_vec(&cell_idx)) - glm::vec2(0.5, 0.5);

            let weights = [
                0.5 * glm::pow(&(glm::vec2(0.5, 0.5) - &cell_diff), &glm::vec2(2.0, 2.0)),
                glm::vec2(0.75, 0.75) - glm::pow(&cell_diff, &glm::vec2(2.0, 2.0)),
                0.5 * glm::pow(&(glm::vec2(0.5, 0.5) + &cell_diff), &glm::vec2(2.0, 2.0)),
            ];

            let mut B: glm::Mat2 = glm::zero();
            for gx in 0..3 {
                for gy in 0..3 {
                    let weight = weights[gx].x * weights[gy].y;
                    let cell_x = glm::vec2(cell_idx.x + gx as u32 - 1,
                                           cell_idx.y + gy as u32 - 1);
                    let cell_idx = cell_x.x as usize * CELL_DIM + cell_x.y as usize;
                    if cell_idx >= self.cells.len() {
                        continue;
                    }

                    let cell_dist = (to_f32_vec(&cell_x) - &p.x) + glm::vec2(0.5, 0.5);
                    let weighted_v = &self.cells[cell_idx].v * weight;

                    let wx = &weighted_v * cell_dist.x;
                    let wy = &weighted_v * cell_dist.y;

                    // TODO Transpose?
                    let term = glm::mat2x2(wx.x, wy.x, wx.y, wy.y);
                    B += term;

                    p.v += weighted_v;
                }
            }

            p.C = B * 4.0;

            p.x += &p.v * dt;
            p.x = glm::clamp(&p.x, 1.0, CELL_DIM as f32-2.0);

            let mut F_new = glm::mat2(1.0, 0.0, 0.0, 1.0);
            F_new += dt * &p.C;
            p.F = &F_new * &p.F;
        }
    }

}