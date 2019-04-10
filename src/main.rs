// struct Vec2D<T> {
//     n_rows: usize,
//     n_cols: usize,
//     data: Vec<T>
// }

// impl<T> Vec2D<T> {
//     fn get(&self, row: usize, col: usize) -> &T {
//          assert!(row < self.n_rows);
//          assert!(col < self.n_cols);
//          &self.data[row * self.n_cols + col]
//     }
// }

struct FluidSquare {
    size: u32,
    dt: f32,
    diff: f32,
    visc: f32,

    s: Vec<f32>,
    density: Vec<f32>,

    vx: Vec<f32>,
    vy: Vec<f32>,

    vx0: Vec<f32>,
    vy0: Vec<f32>,
}

fn index(size: u32, x: u32, y: u32) -> usize {
    (x + y * size) as usize
}

impl FluidSquare {
    fn new(size: u32, diffusion: f32, viscosity: f32, dt: f32) -> FluidSquare {
        let n = size;
        let n_squared = (n * n) as usize;
        FluidSquare {
            size,
            dt,
            diff: diffusion,
            visc: viscosity,
            s: vec![0.0; n_squared],
            density: vec![0.0; n_squared],
            vx: vec![0.0; n_squared],
            vy: vec![0.0; n_squared],
            vx0: vec![0.0; n_squared],
            vy0: vec![0.0; n_squared],
        }
    }

    // fn add_density(&mut self, x: u32, y: u32, amount: f32) {
    //     self.density[index(self.size, x, y)] += amount;
    // }

    // fn add_velocity(&mut self, x: u32, y: u32, amount_x: f32, amount_y: f32) {
    //     let index = index(self.size, x, y);

    //     self.vx[index] += amount_x;
    //     self.vy[index] += amount_y;
    // }
}

fn set_bnd(b: u32, x: &mut Vec<f32>, n: u32) {
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

fn lin_solve(b: u32, x: &mut Vec<f32>, x0: &Vec<f32>, a: f32, c: f32, iter: u32, n: u32) {
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

fn diffuse(b: u32, x: &mut Vec<f32>, x0: &Vec<f32>, diff: f32, dt: f32, iter: u32, n: u32) {
    let a = dt * diff * ((n as f32) - 2.) * ((n as f32) - 2.);
    lin_solve(b, x, x0, a, 1. + 6. * a, iter, n);
}

fn project(
    veloc_x: &mut Vec<f32>,
    veloc_y: &mut Vec<f32>,
    p: &mut Vec<f32>,
    div: &mut Vec<f32>,
    iter: u32,
    n: u32,
) {
    for j in 1..(n - 1) {
        for i in 1..(n - 1) {
            div[index(n, i, j)] = -0.5
                * (veloc_x[index(n, i + 1, j)] - veloc_x[index(n, i - 1, j)]
                    + veloc_y[index(n, i, j + 1)]
                    - veloc_y[index(n, i, j - 1)])
                / (n as f32);
            p[index(n, i, j)] = 0.;
        }
    }
    set_bnd(0, div, n);
    set_bnd(0, p, n);
    lin_solve(0, p, div, 1., 6., iter, n);

    for j in 1..(n - 1) {
        for i in 1..(n - 1) {
            veloc_x[index(n, i, j)] -=
                0.5 * (p[index(n, i + 1, j)] - p[index(n, i - 1, j)]) * (n as f32);
            veloc_y[index(n, i, j)] -=
                0.5 * (p[index(n, i, j + 1)] - p[index(n, i, j - 1)]) * (n as f32);
        }
    }
    set_bnd(1, veloc_x, n);
    set_bnd(2, veloc_y, n);
}

fn advect(
    b: u32,
    d: &mut Vec<f32>,
    d0: &Vec<f32>,
    veloc_x: &Vec<f32>,
    veloc_y: &Vec<f32>,
    dt: f32,
    n: u32,
) {
    let (mut i0, mut i1, mut j0, mut j1);

    let dtx = dt * (n as f32 - 2.);
    let dty = dt * (n as f32 - 2.);

    let (mut s0, mut s1);
    let (mut t0, mut t1);
    let (mut tmp1, mut tmp2);
    let (mut x, mut y);

    let nfloat = n as f32;

    for j in 1..(n - 1) {
        for i in 1..(n - 1) {
            tmp1 = dtx * veloc_x[index(n, i, j)];
            tmp2 = dty * veloc_y[index(n, i, j)];
            x = (i as f32) - tmp1;
            y = (j as f32) - tmp2;

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

    diffuse(1, &mut cube.vx0, &cube.vx, visc, dt, 4, n);
    diffuse(2, &mut cube.vy0, &cube.vy, visc, dt, 4, n);

    project(
        &mut cube.vx0,
        &mut cube.vy0,
        &mut cube.vx,
        &mut cube.vy,
        4,
        n,
    );

    advect(1, &mut cube.vx, &cube.vx0, &cube.vx0, &cube.vy0, dt, n);
    advect(2, &mut cube.vy, &cube.vy0, &cube.vx0, &cube.vy0, dt, n);

    project(
        &mut cube.vx,
        &mut cube.vy,
        &mut cube.vx0,
        &mut cube.vy0,
        4,
        n,
    );

    diffuse(0, &mut cube.s, &cube.density, diff, dt, 4, n);
    advect(0, &mut cube.density, &cube.s, &cube.vx, &cube.vy, dt, n);
}

fn main() {
    let size = 128;
    // let iter = 16;

    let diffusion = 0.2;
    let viscosity = 0.;
    let dt = 0.0000001;
    let mut fluid = FluidSquare::new(size, diffusion, viscosity, dt);

    loop {
        fluid_step(&mut fluid);
    }
}
