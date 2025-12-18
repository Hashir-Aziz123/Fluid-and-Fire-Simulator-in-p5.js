/**
 * FluidBase.js - Abstract Base Class
 * Contains shared velocity solver and grid utilities for both Fluid and Fire.
 * DO NOT instantiate directly - use Fluid or Fire subclasses.
 */

class FluidBase {
    constructor(n, dt, diffusion, viscosity) {
        if (this.constructor === FluidBase) {
            throw new Error("FluidBase is abstract - use Fluid or Fire instead");
        }

        this.N = n; // Grid resolution (inner cells)
        this.dt = dt;
        this.diffusion = diffusion;
        this.viscosity = viscosity;

        // Total size including 1-cell padding on all sides
        this.size = (n + 2) * (n + 2);

        // === SHARED VELOCITY FIELDS ===
        // Both Fluid and Fire need velocity
        this.vX = new Float32Array(this.size);
        this.vY = new Float32Array(this.size);
        this.vXOld = new Float32Array(this.size);
        this.vYOld = new Float32Array(this.size);

        // Working buffers for projection
        this.pressure = new Float32Array(this.size);
        this.divergence = new Float32Array(this.size);
    }

    /**
     * Map 2D coordinates to 1D array index
     */
    idx(x, y) {
        const nx = Math.max(0, Math.min(x, this.N + 1));
        const ny = Math.max(0, Math.min(y, this.N + 1));
        return nx + (this.N + 2) * ny;
    }

    // ========================================
    // VELOCITY SOLVER (Shared by both)
    // ========================================

    /**
     * Solves the velocity field for one timestep
     * @param {number} iterations - Solver iterations for stability
     */
    solveVelocity(iterations) {
        // Diffuse velocity
        this.diffuse(1, this.vXOld, this.vX, this.viscosity, iterations);
        this.diffuse(2, this.vYOld, this.vY, this.viscosity, iterations);
        
        // Make velocity field mass-conserving (incompressible)
        this.project(this.vXOld, this.vYOld, this.vX, this.vY, iterations);

        // Advect velocity through itself
        this.advect(1, this.vX, this.vXOld, this.vXOld, this.vYOld);
        this.advect(2, this.vY, this.vYOld, this.vXOld, this.vYOld);
        
        // Project again for stability
        this.project(this.vX, this.vY, this.vXOld, this.vYOld, iterations);
    }

    /**
     * Diffusion step - spreads quantity via viscosity/diffusion
     */
    diffuse(b, x, x0, diff, iter) {
        const a = this.dt * diff * this.N * this.N;
        this.solveLinearSystem(b, x, x0, a, 1 + 4 * a, iter);
    }

    /**
     * Advection step - moves quantity along velocity field
     */
    advect(b, d, d0, vX, vY) {
        const dt0 = this.dt * this.N;
        
        for (let j = 1; j <= this.N; j++) {
            for (let i = 1; i <= this.N; i++) {
                // Trace particle backward in time
                let x = i - dt0 * vX[this.idx(i, j)];
                let y = j - dt0 * vY[this.idx(i, j)];

                // Clamp to grid boundaries
                if (x < 0.5) x = 0.5; 
                if (x > this.N + 0.5) x = this.N + 0.5;
                if (y < 0.5) y = 0.5; 
                if (y > this.N + 0.5) y = this.N + 0.5;

                // Bilinear interpolation
                const i0 = Math.floor(x); 
                const i1 = i0 + 1;
                const j0 = Math.floor(y); 
                const j1 = j0 + 1;

                const s1 = x - i0; 
                const s0 = 1 - s1;
                const t1 = y - j0; 
                const t0 = 1 - t1;

                d[this.idx(i, j)] = 
                    s0 * (t0 * d0[this.idx(i0, j0)] + t1 * d0[this.idx(i0, j1)]) +
                    s1 * (t0 * d0[this.idx(i1, j0)] + t1 * d0[this.idx(i1, j1)]);
            }
        }
        this.setBoundary(b, d);
    }

    /**
     * Projection step - makes velocity field divergence-free (incompressible)
     */
    project(vX, vY, p, div, iter) {
        // Calculate divergence
        for (let j = 1; j <= this.N; j++) {
            for (let i = 1; i <= this.N; i++) {
                div[this.idx(i, j)] = -0.5 * (
                    vX[this.idx(i + 1, j)] - vX[this.idx(i - 1, j)] +
                    vY[this.idx(i, j + 1)] - vY[this.idx(i, j - 1)]
                ) / this.N;
                p[this.idx(i, j)] = 0;
            }
        }
        this.setBoundary(0, div);
        this.setBoundary(0, p);

        // Solve for pressure
        this.solveLinearSystem(0, p, div, 1, 4, iter);

        // Subtract pressure gradient from velocity
        for (let j = 1; j <= this.N; j++) {
            for (let i = 1; i <= this.N; i++) {
                vX[this.idx(i, j)] -= 0.5 * this.N * 
                    (p[this.idx(i + 1, j)] - p[this.idx(i - 1, j)]);
                vY[this.idx(i, j)] -= 0.5 * this.N * 
                    (p[this.idx(i, j + 1)] - p[this.idx(i, j - 1)]);
            }
        }
        this.setBoundary(1, vX);
        this.setBoundary(2, vY);
    }

    /**
     * Gauss-Seidel iterative solver
     */
    solveLinearSystem(b, x, x0, a, c, iter) {
        const invC = 1.0 / c;
        
        for (let k = 0; k < iter; k++) {
            for (let j = 1; j <= this.N; j++) {
                for (let i = 1; i <= this.N; i++) {
                    x[this.idx(i, j)] = (
                        x0[this.idx(i, j)] + 
                        a * (
                            x[this.idx(i + 1, j)] + 
                            x[this.idx(i - 1, j)] +
                            x[this.idx(i, j + 1)] + 
                            x[this.idx(i, j - 1)]
                        )
                    ) * invC;
                }
            }
            this.setBoundary(b, x);
        }
    }

    /**
     * Boundary conditions
     * b: 0 = continuity, 1 = reflect X, 2 = reflect Y
     */
    setBoundary(b, x) {
        // Edges
        for (let i = 1; i <= this.N; i++) {
            x[this.idx(0, i)] = b === 1 ? -x[this.idx(1, i)] : x[this.idx(1, i)];
            x[this.idx(this.N + 1, i)] = b === 1 ? -x[this.idx(this.N, i)] : x[this.idx(this.N, i)];
            x[this.idx(i, 0)] = b === 2 ? -x[this.idx(i, 1)] : x[this.idx(i, 1)];
            x[this.idx(i, this.N + 1)] = b === 2 ? -x[this.idx(i, this.N)] : x[this.idx(i, this.N)];
        }

        // Corners (average of adjacent edges)
        x[this.idx(0, 0)] = 0.5 * (x[this.idx(1, 0)] + x[this.idx(0, 1)]);
        x[this.idx(0, this.N + 1)] = 0.5 * (x[this.idx(1, this.N + 1)] + x[this.idx(0, this.N)]);
        x[this.idx(this.N + 1, 0)] = 0.5 * (x[this.idx(this.N, 0)] + x[this.idx(this.N + 1, 1)]);
        x[this.idx(this.N + 1, this.N + 1)] = 0.5 * (x[this.idx(this.N, this.N + 1)] + x[this.idx(this.N + 1, this.N)]);
    }

    /**
     * Add velocity at a point (shared interaction)
     */
    addVelocity(x, y, vx, vy) {
        const i = this.idx(x, y);
        this.vX[i] += vx;
        this.vY[i] += vy;
    }

    /**
     * Abstract method - must be implemented by subclasses
     */
    step(params) {
        throw new Error("step() must be implemented by subclass");
    }
}