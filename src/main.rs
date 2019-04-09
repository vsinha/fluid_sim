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

impl FluidSquare {
    fn new(size: u32, diffusion: i32, viscosity: i32, dt: f32) -> FluidSquare {
        let n = size;
        let n_squared = (n * n) as usize;
        FluidSquare {
            size,
            dt,
            diff: diffusion as f32,
            visc: viscosity as f32,
            s: vec![0.0; n_squared],
            density: vec![0.0; n_squared],
            vx: vec![0.0; n_squared],
            vy: vec![0.0; n_squared],
            vx0: vec![0.0; n_squared],
            vy0: vec![0.0; n_squared],
        }
    }
}

fn main() {
    println!("Hello, world!");
}
