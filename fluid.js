class Fluid {
    constructor(dt, diffusion, viscosity) {
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;
        
        // We use the global N from sketch.js, or you can pass it in.
        // For this example, we assume N is accessible or passed. 
        // Let's rely on the global N defined in sketch.js for simplicity 
        // to keep the code matching Day 2's logic.
        this.size = (N + 2) * (N + 2);

        // The "Now" arrays
        this.Vx = new Float32Array(this.size); // Velocity X
        this.Vy = new Float32Array(this.size); // Velocity Y
        this.density = new Float32Array(this.size); // Smoke
        
        // The "Previous" arrays (Buffer)
        this.Vx0 = new Float32Array(this.size);
        this.Vy0 = new Float32Array(this.size);
        this.s = new Float32Array(this.size); // Previous density
    }

    // --- PUBLIC METHODS ---

    step() {
        let visc     = this.visc;
        let diff     = this.diff;
        let dt       = this.dt;
        let Vx      = this.Vx;
        let Vy      = this.Vy;
        let Vx0     = this.Vx0;
        let Vy0     = this.Vy0;
        let s       = this.s;
        let density = this.density;

        // Phase 1: Velocity step
        diffuse(1, Vx0, Vx, visc, dt); // Wind spreads out
        diffuse(2, Vy0, Vy, visc, dt);

        project(Vx0, Vy0, Vx, Vy);  // Fix the wind (swirl)

        advect(1, Vx, Vx0, Vx0, Vy0, dt); // Move the wind via wind
        advect(2, Vy, Vy0, Vx0, Vy0, dt);

        project(Vx, Vy, Vx0, Vy0); // Fix the wind again

        // Phase 2: Density step

        diffuse(0, s, density, diff, dt);  // Smoke spreads out
        advect(0, density, s, Vx, Vy, dt); // Move smoke via wind
    }

    addDensity(x, y, amount) {
        let index = IX(x, y);
        this.density[index] += amount;
    }

    addVelocity(x, y, amountX, amountY) {
        let index = IX(x, y);
        this.Vx[index] += amountX;
        this.Vy[index] += amountY;
    }
}

// --- SOLVER FUNCTIONS (Internal Helpers) ---
// These are outside the class to keep them looking like the original paper's C code,
// but they act on the arrays passed to them.

function IX(x, y) {
    x = constrain(x, 0, N + 1);
    y = constrain(y, 0, N + 1);
    return x + (N + 2) * y;
}

function advect(b, d, d0, velocX, velocY, dt) { // predict cell density via backward trace of flow
    let i0, i1, j0, j1;

    let dtx = dt * (N - 2);
    let dty = dt * (N - 2);

    let s0, s1, t0, t1;
    let tmp1, tmp2, x, y;

    let Nfloat = N;
    let ifloat, jfloat;
    let i, j;

    for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
        for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            tmp1 = dtx * velocX[IX(i, j)];
            tmp2 = dty * velocY[IX(i, j)];
            x = ifloat - tmp1;
            y = jfloat - tmp2;

            if (x < 0.5) x = 0.5;
            if (x > Nfloat + 0.5) x = Nfloat + 0.5;
            i0 = Math.floor(x);
            i1 = i0 + 1.0;
            if (y < 0.5) y = 0.5;
            if (y > Nfloat + 0.5) y = Nfloat + 0.5;
            j0 = Math.floor(y);
            j1 = j0 + 1.0;

            s1 = x - i0;
            s0 = 1.0 - s1;
            t1 = y - j0;
            t0 = 1.0 - t1;

            let i0i = parseInt(i0);
            let i1i = parseInt(i1);
            let j0i = parseInt(j0);
            let j1i = parseInt(j1);

            d[IX(i, j)] =
                s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
                s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
        }
    }
    set_bnd(b, d);
}

function project(velocX, velocY, p, div) { // enforces incompressibliilty for swirl
    for (let j = 1; j < N - 1; j++) {
        for (let i = 1; i < N - 1; i++) {
            div[IX(i, j)] = -0.5 * (
                velocX[IX(i + 1, j)] - velocX[IX(i - 1, j)] +
                velocY[IX(i, j + 1)] - velocY[IX(i, j - 1)]
            ) / N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);

    for (let j = 1; j < N - 1; j++) {
        for (let i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
            velocY[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
        }
    }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

function diffuse(b, x, x0, diff, dt) {
    let a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
}

function lin_solve(b, x, x0, a, c) { // the value of a cell is the average of it's 4 neighbours
    let cRecip = 1.0 / c;
    for (let k = 0; k < iter; k++) {
        for (let j = 1; j < N - 1; j++) {
            for (let i = 1; i < N - 1; i++) {
                x[IX(i, j)] =
                    (x0[IX(i, j)] +
                        a *
                            (x[IX(i + 1, j)] +
                                x[IX(i - 1, j)] +
                                x[IX(i, j + 1)] +
                                x[IX(i, j - 1)])) *
                    cRecip;
            }
        }
        set_bnd(b, x);
    }
}

function set_bnd(b, x) { // keeps the fluid inside the box
    for (let i = 1; i < N - 1; i++) {
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
    }
    for (let j = 1; j < N - 1; j++) {
        x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
        x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
    }

    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N - 1)] = 0.5 * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);
    x[IX(N - 1, 0)] = 0.5 * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);
    x[IX(N - 1, N - 1)] = 0.5 * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);
}

