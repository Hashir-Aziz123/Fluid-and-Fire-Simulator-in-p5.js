/**
 * Fluid simulation using Jos Stam's stable fluids algorithm.
 * Simulates incompressible fluid flow with density, temperature, and velocity fields.
 * Includes Vorticity Confinement for realistic turbulent fire motion.
 */
class Fluid {
    constructor(timeStep, diffusionRate, viscosity) {
        this.timeStep = timeStep;
        this.diffusion = diffusionRate;
        this.viscosity = viscosity;
        
        // Physical properties controlled by UI
        this.buoyancy = 0.0;
        this.coolingRate = 0.99;
        this.velocityDamping = 0.99;
        this.densityFade = 0.995;
        this.turbulenceStrength = 0.0; // Added for Vorticity Confinement strength
        
        // GRID_SIZE + 2 for the 1-cell padding around the simulation grid
        this.gridSize = (GRID_SIZE + 2) * (GRID_SIZE + 2); 
        
        // Density field (visual representation: smoke/dye/ink)
        this.previousDensity = new Float32Array(this.gridSize);
        this.density = new Float32Array(this.gridSize);
        
        // Temperature field (drives buoyancy forces)
        this.previousTemperature = new Float32Array(this.gridSize);
        this.temperature = new Float32Array(this.gridSize);
        
        // Velocity field (fluid motion)
        this.horizontalVelocity = new Float32Array(this.gridSize);
        this.verticalVelocity = new Float32Array(this.gridSize);
        this.previousHorizontalVelocity = new Float32Array(this.gridSize);
        this.previousVerticalVelocity = new Float32Array(this.gridSize);

        // Vorticity fields for turbulence calculation
        this.vorticity = new Float32Array(this.gridSize); 
        this.curlForceX = new Float32Array(this.gridSize);
        this.curlForceY = new Float32Array(this.gridSize);
    }
    
    /**
     * Advances the simulation by one time step.
     */
    step() {
        // --- 1. Forces and Velocity Diffusion ---
        
        // Vprev = V
        this.previousHorizontalVelocity.set(this.horizontalVelocity);
        this.previousVerticalVelocity.set(this.verticalVelocity);

        // Apply Buoyancy forces (hot air rises)
        if (this.buoyancy !== 0) {
            this.applyBuoyancyForce();
        }
        
        // Apply Vorticity Confinement (re-inject turbulence lost to diffusion/advection)
        if (this.turbulenceStrength > 0) {
            this.applyVorticityConfinement();
        }
        
        // Diffuse velocity
        this.diffuseField(1, this.horizontalVelocity, this.previousHorizontalVelocity, this.viscosity);
        this.diffuseField(2, this.verticalVelocity, this.previousVerticalVelocity, this.viscosity);
        
        // --- 2. Projection 1 (Remove divergence from V_diffused) ---
        
        // Vprev = V
        this.previousHorizontalVelocity.set(this.horizontalVelocity);
        this.previousVerticalVelocity.set(this.verticalVelocity);
        
        this.enforceIncompressibility(
            this.horizontalVelocity, 
            this.verticalVelocity,
            this.previousHorizontalVelocity, // Pressure storage
            this.previousVerticalVelocity    // Divergence storage
        );
        
        // --- 3. Advection (Move velocity along itself) ---
        
        // Vprev = V_projected_1
        this.previousHorizontalVelocity.set(this.horizontalVelocity);
        this.previousVerticalVelocity.set(this.verticalVelocity);

        this.advectField(1, this.horizontalVelocity, this.previousHorizontalVelocity,
                         this.previousHorizontalVelocity, this.previousVerticalVelocity);
        this.advectField(2, this.verticalVelocity, this.previousVerticalVelocity,
                         this.previousHorizontalVelocity, this.previousVerticalVelocity);
        
        // --- 4. Projection 2 (Remove divergence from V_advected) ---
        
        this.enforceIncompressibility(
            this.horizontalVelocity, 
            this.verticalVelocity,
            this.previousHorizontalVelocity, // Pressure storage
            this.previousVerticalVelocity    // Divergence storage
        );
        
        // --- 5. Scalar Fields (Temperature and Density) ---
        
        // Temp: Diffuse & Advect
        this.previousTemperature.set(this.temperature);
        this.diffuseField(0, this.temperature, this.previousTemperature, this.diffusion);
        this.advectField(0, this.temperature, this.previousTemperature, 
                         this.horizontalVelocity, this.verticalVelocity);

        // Density: Diffuse & Advect
        this.previousDensity.set(this.density);
        this.diffuseField(0, this.density, this.previousDensity, this.diffusion);
        this.advectField(0, this.density, this.previousDensity,
                         this.horizontalVelocity, this.verticalVelocity);
        
        // --- 6. Decay ---
        this.applyDecay();
    }
    
    /**
     * Adds density (visual smoke/dye) at a specific grid location.
     */
    addDensity(x, y, amount) {
        const index = getGridIndex(x, y);
        this.density[index] += amount;
    }
    
    /**
     * Adds temperature (heat) at a specific grid location.
     */
    addTemperature(x, y, amount) {
        const index = getGridIndex(x, y);
        this.temperature[index] += amount;
    }
    
    /**
     * Adds velocity (motion) at a specific grid location.
     */
    addVelocity(x, y, horizontalAmount, verticalAmount) {
        const index = getGridIndex(x, y);
        this.horizontalVelocity[index] += horizontalAmount;
        this.verticalVelocity[index] += verticalAmount;
    }
    
    /**
     * Applies upward force proportional to temperature AND density (simplified Boussinesq).
     */
    applyBuoyancyForce() {
        // Density factor (for smoke/fuel lift, typically negative)
        const alpha = 0.0001; 
        // Temperature factor (for heat lift, typically positive)
        const beta = this.buoyancy; 
        
        // Loop over the inner, simulated grid cells
        for (let y = 1; y < GRID_SIZE + 1; y++) {
            for (let x = 1; x < GRID_SIZE + 1; x++) {
                const i = getGridIndex(x, y);
                
                // Buoyancy force = (Temp * beta) - (Density * alpha)
                // Negative sign is needed because screen Y increases downward
                const force = -(this.temperature[i] * beta - this.density[i] * alpha);
                
                this.verticalVelocity[i] += force * this.timeStep;
            }
        }
    }
    
    /**
     * Calculates the curl (vorticity) and applies a perpendicular force 
     * to re-inject turbulence, countering numerical damping.
     */
    applyVorticityConfinement() {
        const curl = this.vorticity;
        const normX = this.curlForceX; // Used to store normalized curl gradient (N vector)
        const normY = this.curlForceY; // Used to store normalized curl gradient (N vector)
        
        // Grid cell size in world space (assuming a unit size of 1/GRID_SIZE)
        const h = 1.0 / GRID_SIZE; 
        
        // 1. Calculate Curl (Vorticity) w = dv/dx - du/dy
        for (let y = 1; y < GRID_SIZE + 1; y++) {
            for (let x = 1; x < GRID_SIZE + 1; x++) {
                const i = getGridIndex(x, y);
                
                // Central differencing for derivatives
                const du_dy = (this.horizontalVelocity[getGridIndex(x, y + 1)] - this.horizontalVelocity[getGridIndex(x, y - 1)]) / (2 * h);
                const dv_dx = (this.verticalVelocity[getGridIndex(x + 1, y)] - this.verticalVelocity[getGridIndex(x - 1, y)]) / (2 * h);
                
                curl[i] = dv_dx - du_dy;
            }
        }
        
        // 2. Calculate the "N" vector (Normalized Gradient of Curl Magnitude) N = grad(|w|) / |grad(|w|)|
        for (let y = 1; y < GRID_SIZE + 1; y++) {
            for (let x = 1; x < GRID_SIZE + 1; x++) {
                const i = getGridIndex(x, y);
                
                // Central differencing for the gradient of |curl|
                const abs_curl_x = (Math.abs(curl[getGridIndex(x + 1, y)]) - Math.abs(curl[getGridIndex(x - 1, y)])) / (2 * h);
                const abs_curl_y = (Math.abs(curl[getGridIndex(x, y + 1)]) - Math.abs(curl[getGridIndex(x, y - 1)])) / (2 * h);
                
                // Magnitude of the gradient
                let magnitude = Math.sqrt(abs_curl_x * abs_curl_x + abs_curl_y * abs_curl_y) + 1e-30; 
                
                // Normalize N
                normX[i] = abs_curl_x / magnitude;
                normY[i] = abs_curl_y / magnitude;
            }
        }
        
        // 3. Apply Confinement Force (F = N x w)
        for (let i = 0; i < this.gridSize; i++) {
            // Force is perpendicular to N and proportional to curl magnitude
            // Fx = Ny * w * strength
            const forceX = normY[i] * curl[i] * this.turbulenceStrength; 
            // Fy = -Nx * w * strength
            const forceY = -normX[i] * curl[i] * this.turbulenceStrength; 
            
            // Add force to velocity field
            this.horizontalVelocity[i] += forceX * this.timeStep;
            this.verticalVelocity[i] += forceY * this.timeStep;
        }
        
        this.setBoundaryConditions(1, this.horizontalVelocity);
        this.setBoundaryConditions(2, this.verticalVelocity);
    }
    
    /**
     * Applies natural decay to all fields over time.
     * Velocity dampens, density fades, temperature cools.
     */
    applyDecay() {
        for (let i = 0; i < this.gridSize; i++) {
            // Velocity friction (dampens flow)
            this.horizontalVelocity[i] *= this.velocityDamping;
            this.verticalVelocity[i] *= this.velocityDamping;
            
            // Density fade (prevents smoke/dye from filling the screen)
            this.density[i] *= this.densityFade;
            
            // Temperature cooling (critical for defining the flame shape)
            if (this.buoyancy > 0) {
                this.temperature[i] *= this.coolingRate;
            } else {
                this.temperature[i] = 0;
            }
        }
    }
    
    // --- Jos Stam's Stable Fluids Core Methods ---

    /**
     * Diffuses a field (spreads values to neighbors). Solved iteratively.
     */
    diffuseField(boundaryType, target, source, diffusionRate) {
        // Pre-calculate factor A based on the standard Stable Fluids formulation
        const diffusionFactor = this.timeStep * diffusionRate * (GRID_SIZE) * (GRID_SIZE); 
        // Factor C is used for the Gauss-Seidel relaxation (1 + 6 * A)
        this.solveLinearSystem(boundaryType, target, source, diffusionFactor, 1 + 6 * diffusionFactor);
    }
    
    /**
     * Advects a field (moves values along velocity field using semi-Lagrangian).
     */
    advectField(boundaryType, target, source, velocityX, velocityY) {
        // Time scaling factor for the advection distance
        const timeScale = this.timeStep * GRID_SIZE; 
        
        for (let y = 1; y < GRID_SIZE + 1; y++) {
            for (let x = 1; x < GRID_SIZE + 1; x++) {
                const index = getGridIndex(x, y);
                
                // Trace particle backward in time: p(t-dt) = p(t) - V(p(t)) * dt
                let sourceX = x - timeScale * velocityX[index];
                let sourceY = y - timeScale * velocityY[index];
                
                // Clamp to valid grid range (excluding boundary cells for interpolation)
                sourceX = constrain(sourceX, 0.5, GRID_SIZE + 0.5);
                sourceY = constrain(sourceY, 0.5, GRID_SIZE + 0.5);
                
                // Bilinear interpolation indices
                const x0 = Math.floor(sourceX);
                const x1 = x0 + 1;
                const y0 = Math.floor(sourceY);
                const y1 = y0 + 1;
                
                // Interpolation weights (x axis)
                const weightX1 = sourceX - x0;
                const weightX0 = 1.0 - weightX1;
                // Interpolation weights (y axis)
                const weightY1 = sourceY - y0;
                const weightY0 = 1.0 - weightY1;
                
                // Bilinear Interpolation calculation
                target[index] = 
                    weightX0 * (weightY0 * source[getGridIndex(x0, y0)] + 
                                weightY1 * source[getGridIndex(x0, y1)]) +
                    weightX1 * (weightY0 * source[getGridIndex(x1, y0)] + 
                                weightY1 * source[getGridIndex(x1, y1)]);
            }
        }
        
        this.setBoundaryConditions(boundaryType, target);
    }
    
    /**
     * Projects velocity field to be divergence-free (incompressible) using pressure.
     * Solved iteratively (part of the core stable fluids).
     */
    enforceIncompressibility(velocityX, velocityY, pressure, divergence) {
        // 1. Compute Divergence (div V = partial_u/partial_x + partial_v/partial_y)
        for (let y = 1; y < GRID_SIZE + 1; y++) {
            for (let x = 1; x < GRID_SIZE + 1; x++) {
                const index = getGridIndex(x, y);
                
                // Divergence (D = -h * div V) using central differencing
                divergence[index] = -0.5 * (
                    velocityX[getGridIndex(x + 1, y)] - velocityX[getGridIndex(x - 1, y)] +
                    velocityY[getGridIndex(x, y + 1)] - velocityY[getGridIndex(x, y - 1)]
                ) * GRID_SIZE; // Multiply by grid size to account for cell spacing
                
                pressure[index] = 0; // Initialize pressure to zero
            }
        }
        
        this.setBoundaryConditions(0, divergence);
        this.setBoundaryConditions(0, pressure);
        
        // 2. Solve Poisson Equation for Pressure (p = div V)
        this.solveLinearSystem(0, pressure, divergence, 1, 6); // 1 and 6 are standard factors
        
        // 3. Subtract Pressure Gradient from Velocity (V = V - grad p)
        for (let y = 1; y < GRID_SIZE + 1; y++) {
            for (let x = 1; x < GRID_SIZE + 1; x++) {
                const index = getGridIndex(x, y);
                
                // Subtract pressure gradient (grad P) using central differencing
                velocityX[index] -= 0.5 * (
                    pressure[getGridIndex(x + 1, y)] - pressure[getGridIndex(x - 1, y)]
                ) * GRID_SIZE;
                
                velocityY[index] -= 0.5 * (
                    pressure[getGridIndex(x, y + 1)] - pressure[getGridIndex(x, y - 1)]
                ) * GRID_SIZE;
            }
        }
        
        this.setBoundaryConditions(1, velocityX);
        this.setBoundaryConditions(2, velocityY);
    }
    
    /**
     * Solves a linear system Ax = b using Gauss-Seidel relaxation (iterative solver).
     */
    solveLinearSystem(boundaryType, target, source, factorA, factorC) {
        const inverseFactor = 1.0 / factorC;
        const iterations = SOLVER_ITERATIONS || 4; // Use global const SOLVER_ITERATIONS
        
        for (let iteration = 0; iteration < iterations; iteration++) {
            for (let y = 1; y < GRID_SIZE + 1; y++) {
                for (let x = 1; x < GRID_SIZE + 1; x++) {
                    const index = getGridIndex(x, y);
                    
                    // Gauss-Seidel update: x_ij = (b_ij + A * sum(x_neighbors)) / C
                    target[index] = (
                        source[index] +
                        factorA * (
                            target[getGridIndex(x + 1, y)] +
                            target[getGridIndex(x - 1, y)] +
                            target[getGridIndex(x, y + 1)] +
                            target[getGridIndex(x, y - 1)]
                        )
                    ) * inverseFactor;
                }
            }
            
            this.setBoundaryConditions(boundaryType, target);
        }
    }
    
    /**
     * Sets boundary conditions for the field (Handles no-slip/no-through boundaries).
     */
    setBoundaryConditions(boundaryType, field) {
        // Helper to determine sign for boundary conditions
        const signX = boundaryType === 1 ? -1 : 1;
        const signY = boundaryType === 2 ? -1 : 1;
        
        // Top and bottom edges
        for (let x = 1; x < GRID_SIZE + 1; x++) {
            field[getGridIndex(x, 0)] = signY * field[getGridIndex(x, 1)];
            field[getGridIndex(x, GRID_SIZE + 1)] = signY * field[getGridIndex(x, GRID_SIZE)];
        }
        
        // Left and right edges
        for (let y = 1; y < GRID_SIZE + 1; y++) {
            field[getGridIndex(0, y)] = signX * field[getGridIndex(1, y)];
            field[getGridIndex(GRID_SIZE + 1, y)] = signX * field[getGridIndex(GRID_SIZE, y)];
        }
        
        // Corners (average of adjacent cells)
        field[getGridIndex(0, 0)] = 0.5 * (field[getGridIndex(1, 0)] + field[getGridIndex(0, 1)]);
        field[getGridIndex(0, GRID_SIZE + 1)] = 0.5 * (field[getGridIndex(1, GRID_SIZE + 1)] + field[getGridIndex(0, GRID_SIZE)]);
        field[getGridIndex(GRID_SIZE + 1, 0)] = 0.5 * (field[getGridIndex(GRID_SIZE, 0)] + field[getGridIndex(GRID_SIZE + 1, 1)]);
        field[getGridIndex(GRID_SIZE + 1, GRID_SIZE + 1)] = 0.5 * (field[getGridIndex(GRID_SIZE, GRID_SIZE + 1)] + field[getGridIndex(GRID_SIZE + 1, GRID_SIZE)]);
    }
}

// Helper function for grid indexing (needed globally for sketch.js)
function getGridIndex(x, y) {
    x = constrain(x, 0, GRID_SIZE + 1);
    y = constrain(y, 0, GRID_SIZE + 1);
    return x + (GRID_SIZE + 2) * y;
}