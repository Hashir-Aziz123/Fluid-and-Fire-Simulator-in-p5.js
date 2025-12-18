/**
 * Fire.js - Combustion Simulation
 * Updated with Vorticity Confinement for "Swirl"
 */

class Fire extends FluidBase {
    constructor(n, dt, diffusion, viscosity) {
        super(n, dt, diffusion, viscosity);

        // === FIRE-SPECIFIC FIELDS ===
        this.fuel = new Float32Array(this.size);
        this.fuelOld = new Float32Array(this.size);
        
        this.temperature = new Float32Array(this.size);
        this.tempOld = new Float32Array(this.size);
        
        this.smoke = new Float32Array(this.size);
        this.smokeOld = new Float32Array(this.size);
        
        this.reaction = new Float32Array(this.size);

        // === NEW: Vorticity Buffer (The Swirl Storage) ===
        this.curl = new Float32Array(this.size);

        // === FIRE PHYSICS CONSTANTS ===
        this.ignitionTemp = 20;        
        this.burnRate = 0.003;         
        this.heatRelease = 15;         
        this.smokeGeneration = 2.5;    
        this.heatDiffusion = 0.02;     
        this.ambientCooling = 0.98;    
    }

    /**
     * Main simulation step for fire
     */
    step(params) {
        // 1. COMBUSTION REACTION 
        this.combustion(params);

        // 2. HEAT DIFFUSION 
        this.diffuseHeat(params.iterations);

        // 3. TEMPERATURE-DRIVEN BUOYANCY
        for (let i = 0; i < this.size; i++) {
            // Hot air rises
            this.vY[i] -= this.temperature[i] * params.buoyancy * 0.0003;
        }

        // === NEW: APPLY SWIRL (Vorticity) ===
        // We use a high default strength if not provided to force the curl
        if (params.vorticity > 0) {
            this.applyVorticityConfinement(params.vorticity);
        }

        // 4. SOLVE VELOCITY FIELD
        this.solveVelocity(params.iterations);

        // 5. ADVECT SCALARS
        this.diffuse(0, this.fuelOld, this.fuel, this.diffusion * 0.1, params.iterations);
        this.advect(0, this.fuel, this.fuelOld, this.vX, this.vY);

        this.advect(0, this.temperature, this.tempOld, this.vX, this.vY);

        this.diffuse(0, this.smokeOld, this.smoke, this.diffusion, params.iterations);
        this.advect(0, this.smoke, this.smokeOld, this.vX, this.vY);

        // 6. COOLING AND DECAY
        for (let i = 0; i < this.size; i++) {
            this.vX[i] *= params.damping;
            this.vY[i] *= params.damping;
            this.temperature[i] *= params.cooling;
            this.smoke[i] *= params.fade;
        }
    }

    /**
     * === NEW FUNCTION: Calculates and applies spin ===
     */
    applyVorticityConfinement(strength) {
        // Calculate curl (vorticity)
        for (let j = 1; j <= this.N; j++) {
            for (let i = 1; i <= this.N; i++) {
                const dvx_dy = (this.vX[this.idx(i, j + 1)] - this.vX[this.idx(i, j - 1)]) * 0.5;
                const dvy_dx = (this.vY[this.idx(i + 1, j)] - this.vY[this.idx(i - 1, j)]) * 0.5;
                this.curl[this.idx(i, j)] = dvy_dx - dvx_dy;
            }
        }

        // Apply force to increase swirl
        for (let j = 2; j < this.N; j++) {
            for (let i = 2; i < this.N; i++) {
                const dx = (Math.abs(this.curl[this.idx(i + 1, j)]) - Math.abs(this.curl[this.idx(i - 1, j)])) * 0.5;
                const dy = (Math.abs(this.curl[this.idx(i, j + 1)]) - Math.abs(this.curl[this.idx(i, j - 1)])) * 0.5;

                const length = Math.sqrt(dx * dx + dy * dy) + 1e-5;
                const nx = dx / length;
                const ny = dy / length;

                const curlValue = this.curl[this.idx(i, j)];

                this.vX[this.idx(i, j)] += ny * curlValue * strength * this.dt;
                this.vY[this.idx(i, j)] -= nx * curlValue * strength * this.dt;
            }
        }
    }

    // ... (Keep combustion, diffuseHeat, addFuel, addTemperature, getIntensity same as before) ...
    
    combustion(params) {
       // Paste your original combustion code here (no changes needed)
       for (let j = 1; j <= this.N; j++) {
            for (let i = 1; i <= this.N; i++) {
                const idx = this.idx(i, j);
                const f = this.fuel[idx];
                const t = this.temperature[idx];
                if (f > 1.0 && t > this.ignitionTemp) {
                    const fuelFactor = Math.min(f / 100, 1.0);
                    const tempFactor = Math.min((t - this.ignitionTemp) / 200, 1.0);
                    const reactionRate = fuelFactor * tempFactor;
                    this.reaction[idx] = reactionRate;
                    this.fuel[idx] -= reactionRate * this.burnRate;
                    if (this.fuel[idx] < 0) this.fuel[idx] = 0;
                    this.temperature[idx] += reactionRate * this.heatRelease;
                    this.smoke[idx] += reactionRate * this.smokeGeneration;
                } else {
                    this.reaction[idx] = 0;
                }
            }
        }
    }

    diffuseHeat(iterations) {
        // Paste your original diffuseHeat code here (no changes needed)
        for (let i = 0; i < this.size; i++) {
            this.tempOld[i] = this.temperature[i];
        }
        const a = this.dt * this.heatDiffusion * this.N * this.N;
        this.solveLinearSystem(0, this.temperature, this.tempOld, a, 1 + 4 * a, iterations);
    }

    addFuel(x, y, amount) { this.fuel[this.idx(x, y)] += amount; }
    addTemperature(x, y, amount) { this.temperature[this.idx(x, y)] += amount; }
    getIntensity(x, y) {
        const idx = this.idx(x, y);
        return { temperature: this.temperature[idx], smoke: this.smoke[idx], fuel: this.fuel[idx], reaction: this.reaction[idx] };
    }
}