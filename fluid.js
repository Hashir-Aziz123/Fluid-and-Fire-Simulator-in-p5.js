/**
 * Fluid.js - Water/Smoke Simulation
 * Extends FluidBase with density field and vorticity effects.
 */

class Fluid extends FluidBase {
    constructor(n, dt, diffusion, viscosity) {
        super(n, dt, diffusion, viscosity);

        // === FLUID-SPECIFIC FIELDS ===
        this.density = new Float32Array(this.size);
        this.densityOld = new Float32Array(this.size);
        
        // Vorticity buffer for swirl effects
        this.curl = new Float32Array(this.size);
    }

    /**
     * Main simulation step for fluid
     */
    step(params) {
        // 1. Apply vorticity confinement (swirl effect)
        if (params.vorticity > 0) {
            this.applyVorticityConfinement(params.vorticity);
        }

        // 2. Solve velocity field
        this.solveVelocity(params.iterations);

        // 3. Advect and diffuse density
        this.diffuse(0, this.densityOld, this.density, this.diffusion, params.iterations);
        this.advect(0, this.density, this.densityOld, this.vX, this.vY);

        // 4. Apply decay effects
        for (let i = 0; i < this.size; i++) {
            this.vX[i] *= params.damping;
            this.vY[i] *= params.damping;
            this.density[i] *= params.fade;
        }
    }

    /**
     * Vorticity Confinement - adds rotational energy to prevent numerical dissipation
     */
    applyVorticityConfinement(strength) {
        // Calculate curl (vorticity) at each cell
        for (let j = 1; j <= this.N; j++) {
            for (let i = 1; i <= this.N; i++) {
                const dvx_dy = (this.vX[this.idx(i, j + 1)] - this.vX[this.idx(i, j - 1)]) * 0.5;
                const dvy_dx = (this.vY[this.idx(i + 1, j)] - this.vY[this.idx(i - 1, j)]) * 0.5;
                this.curl[this.idx(i, j)] = dvy_dx - dvx_dy;
            }
        }

        // Apply confinement force based on gradient of curl magnitude
        for (let j = 2; j < this.N; j++) {
            for (let i = 2; i < this.N; i++) {
                const dx = (Math.abs(this.curl[this.idx(i + 1, j)]) - 
                           Math.abs(this.curl[this.idx(i - 1, j)])) * 0.5;
                const dy = (Math.abs(this.curl[this.idx(i, j + 1)]) - 
                           Math.abs(this.curl[this.idx(i, j - 1)])) * 0.5;

                const length = Math.sqrt(dx * dx + dy * dy) + 1e-5;
                const nx = dx / length;
                const ny = dy / length;

                const curlValue = this.curl[this.idx(i, j)];

                // Add perpendicular force to velocity
                this.vX[this.idx(i, j)] += ny * curlValue * strength * this.dt;
                this.vY[this.idx(i, j)] -= nx * curlValue * strength * this.dt;
            }
        }
    }

    /**
     * Add density at a point (user interaction)
     */
    addDensity(x, y, amount) {
        this.density[this.idx(x, y)] += amount;
    }
}