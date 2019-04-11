extern crate glutin_window;
extern crate graphics;
extern crate opengl_graphics;
extern crate palette;
extern crate piston;

use glutin_window::GlutinWindow as Window;
use opengl_graphics::{GlGraphics, OpenGL};
use piston::event_loop::*;
use piston::input::*;
use piston::window::WindowSettings;

struct FluidSquare {
    size: u32,
    dt: f64,
    diff: f64,
    visc: f64,
    iter: u32,

    s: Vec<f64>,
    density: Vec<f64>,

    vx: Vec<f64>,
    vy: Vec<f64>,

    vx0: Vec<f64>,
    vy0: Vec<f64>,
}

fn index(size: u32, x: u32, y: u32) -> usize {
    let x = constrain(x, 0, size - 1);
    let y = constrain(y, 0, size - 1);
    (x + y * size) as usize
}

fn constrain<T: PartialOrd>(val: T, min: T, max: T) -> T {
    if val < min {
        min
    } else if val >= max {
        max
    } else {
        val
    }
}

impl FluidSquare {
    fn new(size: u32, diffusion: f64, viscosity: f64, dt: f64, iter: u32) -> FluidSquare {
        let n = size;
        let n_squared = (n * n) as usize;
        FluidSquare {
            size,
            dt,
            diff: diffusion,
            visc: viscosity,
            iter,
            s: vec![0.0; n_squared],
            density: vec![0.0; n_squared],
            vx: vec![0.0; n_squared],
            vy: vec![0.0; n_squared],
            vx0: vec![0.0; n_squared],
            vy0: vec![0.0; n_squared],
        }
    }

    fn add_density(&mut self, x: u32, y: u32, amount: f64) {
        self.density[index(self.size, x, y)] += amount;
    }

    fn add_velocity(&mut self, x: u32, y: u32, amount_x: f64, amount_y: f64) {
        let index = index(self.size, x, y);

        self.vx[index] += amount_x;
        self.vy[index] += amount_y;
    }
}

fn set_bnd(b: u32, x: &mut Vec<f64>, n: u32) {
    for i in 1..(n - 1) {
        if b == 2 {
            x[index(n, i, 0)] = -x[index(n, i, 1)];
            x[index(n, i, n - 1)] = -x[index(n, i, n - 2)];
        } else {
            x[index(n, i, 0)] = x[index(n, i, 1)];
            x[index(n, i, n - 1)] = x[index(n, i, n - 2)];
        }
    }

    for j in 1..(n - 1) {
        if b == 2 {
            x[index(n, 0, j)] = -x[index(n, 1, j)];
            x[index(n, n - 1, j)] = -x[index(n, n - 2, j)];
        } else {
            x[index(n, 0, j)] = x[index(n, 1, j)];
            x[index(n, n - 1, j)] = x[index(n, n - 2, j)];
        }
    }

    x[index(n, 0, 0)] = 0.5 * (x[index(n, 1, 0)] + x[index(n, 0, 1)]);
    x[index(n, 0, n - 1)] = 0.5 * (x[index(n, 1, n - 1)] + x[index(n, 0, n - 2)]);
    x[index(n, n - 1, 0)] = 0.5 * (x[index(n, n - 2, 0)] + x[index(n, n - 1, 1)]);
    x[index(n, n - 1, n - 1)] = 0.5 * (x[index(n, n - 2, n - 1)] + x[index(n, n - 1, n - 2)]);
}

fn lin_solve(b: u32, x: &mut Vec<f64>, x0: &Vec<f64>, a: f64, c: f64, iter: u32, n: u32) {
    let c_recip = 1.0 / c;
    for _k in 0..iter {
        for j in 1..(n - 1) {
            for i in 1..(n - 1) {
                x[index(n, i, j)] = (x0[index(n, i, j)]
                    + a * (x[index(n, i + 1, j)]
                        + x[index(n, i - 1, j)]
                        + x[index(n, i, j + 1)]
                        + x[index(n, i, j - 1)]
                        + x[index(n, i, j + 1)]
                        + x[index(n, i, j - 1)]))
                    * c_recip;
            }
        }
        set_bnd(b, x, n);
    }
}

fn diffuse(b: u32, x: &mut Vec<f64>, x0: &Vec<f64>, diff: f64, dt: f64, iter: u32, n: u32) {
    let a = dt * diff * ((n as f64) - 2.) * ((n as f64) - 2.);
    lin_solve(b, x, x0, a, 1. + 6. * a, iter, n);
}

fn project(
    veloc_x: &mut Vec<f64>,
    veloc_y: &mut Vec<f64>,
    p: &mut Vec<f64>,
    div: &mut Vec<f64>,
    iter: u32,
    n: u32,
) {
    for j in 1..(n - 1) {
        for i in 1..(n - 1) {
            div[index(n, i, j)] = -0.5
                * (veloc_x[index(n, i + 1, j)] - veloc_x[index(n, i - 1, j)]
                    + veloc_y[index(n, i, j + 1)]
                    - veloc_y[index(n, i, j - 1)])
                / (n as f64);
            p[index(n, i, j)] = 0.;
        }
    }
    set_bnd(0, div, n);
    set_bnd(0, p, n);
    lin_solve(0, p, div, 1., 6., iter, n);

    for j in 1..(n - 1) {
        for i in 1..(n - 1) {
            veloc_x[index(n, i, j)] -=
                0.5 * (p[index(n, i + 1, j)] - p[index(n, i - 1, j)]) * (n as f64);
            veloc_y[index(n, i, j)] -=
                0.5 * (p[index(n, i, j + 1)] - p[index(n, i, j - 1)]) * (n as f64);
        }
    }
    set_bnd(1, veloc_x, n);
    set_bnd(2, veloc_y, n);
}

fn advect(
    b: u32,
    d: &mut Vec<f64>,
    d0: &Vec<f64>,
    veloc_x: &Vec<f64>,
    veloc_y: &Vec<f64>,
    dt: f64,
    n: u32,
) {
    let (mut i0, mut i1, mut j0, mut j1);

    let dtx = dt * (n as f64 - 2.);
    let dty = dt * (n as f64 - 2.);

    let (mut s0, mut s1);
    let (mut t0, mut t1);
    let (mut tmp1, mut tmp2);
    let (mut x, mut y);

    let nfloat = n as f64;

    for j in 1..(n - 1) {
        for i in 1..(n - 1) {
            tmp1 = dtx * veloc_x[index(n, i, j)];
            tmp2 = dty * veloc_y[index(n, i, j)];
            x = (i as f64) - tmp1;
            y = (j as f64) - tmp2;

            if x < 0.5 {
                x = 0.5
            };
            if x > (nfloat + 0.5) {
                x = nfloat + 0.5
            };
            i0 = x.floor();
            i1 = i0 + 1.0;
            if y < 0.5 {
                y = 0.5
            };
            if y > (nfloat + 0.5) {
                y = nfloat + 0.5
            };
            j0 = y.floor();
            j1 = j0 + 1.0;

            s1 = x - i0;
            s0 = 1.0 - s1;
            t1 = y - j0;
            t0 = 1.0 - t1;

            let i0i = i0 as u32;
            let i1i = i1 as u32;
            let j0i = j0 as u32;
            let j1i = j1 as u32;

            d[index(n, i, j)] = s0 * t0 * d0[index(n, i0i, j0i)]
                + t1 * d0[index(n, i0i, j1i)]
                + s1 * t0 * d0[index(n, i1i, j0i)]
                + t1 * d0[index(n, i1i, j1i)];
        }
    }
    set_bnd(b, d, n);
}

fn fluid_step(cube: &mut FluidSquare) {
    let n = cube.size;
    let visc = cube.visc;
    let diff = cube.diff;
    let dt = cube.dt;
    let iter = cube.iter;

    diffuse(1, &mut cube.vx0, &cube.vx, visc, dt, iter, n);
    diffuse(2, &mut cube.vy0, &cube.vy, visc, dt, iter, n);

    project(
        &mut cube.vx0,
        &mut cube.vy0,
        &mut cube.vx,
        &mut cube.vy,
        iter,
        n,
    );

    advect(1, &mut cube.vx, &cube.vx0, &cube.vx0, &cube.vy0, dt, n);
    advect(2, &mut cube.vy, &cube.vy0, &cube.vx0, &cube.vy0, dt, n);

    project(
        &mut cube.vx,
        &mut cube.vy,
        &mut cube.vx0,
        &mut cube.vy0,
        iter,
        n,
    );

    diffuse(0, &mut cube.s, &cube.density, diff, dt, iter, n);
    advect(0, &mut cube.density, &cube.s, &cube.vx, &cube.vy, dt, n);
}

pub struct App {
    gl: GlGraphics,
    fluid: FluidSquare,
    width: u32,
    height: u32,
    scale: u32,
}

impl App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        let fluid = &self.fluid;
        let square = rectangle::square(0.0, 0.0, self.scale as f64);
        let scale = self.scale;

        self.gl.draw(args.viewport(), |c, gl| {
            // Clear the screen.
            clear([0., 0., 0., 0.0], gl);

            let n = fluid.size;
            for i in 0..n {
                for j in 0..n {
                    let x = i * scale;
                    let y = j * scale;
                    let transform = c.transform.trans(x as f64, y as f64);

                    let d = fluid.density[index(n, i, j)];
                    let h = (d + 50. % 255.) as f32;
                    let s = 200. as f32;
                    let v = d as f32;
                    let hsv = palette::Hsv::new(h, s, v);
                    let rgb: palette::rgb::Rgb = palette::rgb::Rgb::from(hsv);
                    let color = [rgb.red, rgb.green, rgb.blue, 1.];

                    rectangle(color, square, transform, gl);
                }
            }
        });
    }

    fn update(&mut self, _args: &UpdateArgs) {
        let cx = (0.5 * (self.width / self.scale) as f64) as u32;
        let cy = (0.5 * (self.height / self.scale) as f64) as u32;

        FluidSquare::add_density(&mut self.fluid, cx, cy, 0.00001);
        FluidSquare::add_velocity(&mut self.fluid, cx, cy, 1., 1.);

        fluid_step(&mut self.fluid);
    }
}

fn main() {
    // Change this to OpenGL::V2_1 if not working.
    let opengl = OpenGL::V3_2;

    let size = 128;
    let scale = 4;
    let width = size * scale;
    let height = size * scale;

    // Create an Glutin window.
    let mut window: Window = WindowSettings::new("Something", [width, height])
        .opengl(opengl)
        .exit_on_esc(true)
        .build()
        .unwrap();

    let iter = 4;
    let diffusion = 0.2;
    let viscosity = 0.;
    let dt = 0.000001;

    // Create a new game and run it.
    let mut app = App {
        gl: GlGraphics::new(opengl),
        fluid: FluidSquare::new(size, diffusion, viscosity, dt, iter),
        width,
        height,
        scale,
    };

    let mut events = Events::new(EventSettings::new());
    while let Some(e) = events.next(&mut window) {
        if let Some(r) = e.render_args() {
            app.render(&r);
        }

        if let Some(u) = e.update_args() {
            app.update(&u);
        }
    }
}
