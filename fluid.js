class Fluid {
    constructor(dt, diffusion, viscosity) {
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;
        
        // --- CONFIGURATION PROPERTIES ---
        this.buoyancy = 0.0; 
        this.cooling = 0.99; 
        
        this.size = (N + 2) * (N + 2);

        // 1. DENSITY (The Ink/Smoke)
        this.s = new Float32Array(this.size);
        this.density = new Float32Array(this.size);
        
        // 2. TEMPERATURE (The Heat) <--- NEW
        this.t0 = new Float32Array(this.size);
        this.t = new Float32Array(this.size);
        
        // 3. VELOCITY (The Wind)
        this.Vx = new Float32Array(this.size);
        this.Vy = new Float32Array(this.size);
        this.Vx0 = new Float32Array(this.size);
        this.Vy0 = new Float32Array(this.size);
    }

    step() {
        let visc = this.visc;
        let diff = this.diff;
        let dt   = this.dt;
        let Vx   = this.Vx;
        let Vy   = this.Vy;
        let Vx0  = this.Vx0;
        let Vy0  = this.Vy0;
        let s    = this.s;
        let density = this.density;
        let t    = this.t;    // <--- NEW
        let t0   = this.t0;   // <--- NEW

        // A. SOLVE TEMPERATURE (Heat spreads and moves)
        // Heat diffuses slightly faster than smoke usually, but we use same 'diff' for simplicity
        diffuse(0, t0, t, diff, dt);
        advect(0, t, t0, Vx, Vy, dt);

        // B. APPLY FORCES (Buoyancy based on Temp)
        if (this.buoyancy > 0) {
            this.applyBuoyancy();
        }

        // C. SOLVE VELOCITY
        diffuse(1, Vx0, Vx, visc, dt);
        diffuse(2, Vy0, Vy, visc, dt);
        project(Vx0, Vy0, Vx, Vy);
        advect(1, Vx, Vx0, Vx0, Vy0, dt);
        advect(2, Vy, Vy0, Vx0, Vy0, dt);
        project(Vx, Vy, Vx0, Vy0);

        // D. SOLVE DENSITY
        diffuse(0, s, density, diff, dt);
        advect(0, density, s, Vx, Vy, dt);
        
        // E. COOLING & DECAY
        this.decay();
    }

    // --- HELPER METHODS ---

    addDensity(x, y, amount) {
        let index = IX(x, y);
        this.density[index] += amount;
    }
    
    // NEW: Add Heat
    addTemperature(x, y, amount) {
        let index = IX(x, y);
        this.t[index] += amount;
    }

    addVelocity(x, y, amountX, amountY) {
        let index = IX(x, y);
        this.Vx[index] += amountX;
        this.Vy[index] += amountY;
    }

    // UPDATED: Simulates heat rising based on TEMPERATURE
    applyBuoyancy() {
        for (let i = 0; i < this.size; i++) {
            let temp = this.t[i];
            
            // Only hot air rises
            if (temp > 0) {
                // "Up" is Negative Y. 
                // We use temperature to drive the force.
                this.Vy[i] -= temp * this.buoyancy;
            }
        }
    }

    // UPDATED: Decay (No more hacks!)
    decay() {
        for (let i = 0; i < this.size; i++) {
            // 1. Velocity Damping
            this.Vx[i] *= 0.99; 
            this.Vy[i] *= 0.99;

            // 2. Density Decay (Slow fade)
            this.density[i] *= 0.995; 

            // 3. Temperature Cooling (Fast fade) <--- NEW
            // Heat disappears much faster than smoke!
            if (this.buoyancy > 0) {
                 this.t[i] *= this.cooling; 
            } else {
                 this.t[i] = 0; // Clear heat if not in fire mode
            }
        }
    }
}

// ... (KEEP ALL THE SOLVER FUNCTIONS: IX, advect, project, diffuse, lin_solve, set_bnd) ...
// ... (They are unchanged) ...

// --- STANDARD JOS STAM SOLVER FUNCTIONS (Unchanged) ---
// (These remain exactly the same as previous steps)

function IX(x, y) {
    x = constrain(x, 0, N + 1);
    y = constrain(y, 0, N + 1);
    return x + (N + 2) * y;
}

function advect(b, d, d0, velocX, velocY, dt) {
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

function project(velocX, velocY, p, div) {
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

function lin_solve(b, x, x0, a, c) {
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

function set_bnd(b, x) {
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
