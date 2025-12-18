/**
 * Renderer.js - Visual Artist
 * Separate rendering paths for Fluid and Fire simulations.
 */

class SimulationRenderer {
    constructor(p5Context, n, scale) {
        this.p = p5Context;
        this.N = n;
        this.scale = scale;
        this.buffer = this.p.createImage(n, n);
    }

    /**
     * Render fluid simulation (density-based)
     */
    renderFluid(fluid, params) {
        this.buffer.loadPixels();

        for (let y = 1; y <= this.N; y++) {
            for (let x = 1; x <= this.N; x++) {
                const i = fluid.idx(x, y);
                const density = fluid.density[i];
                
                const pxIdx = 4 * ((x - 1) + (y - 1) * this.N);

                // Simple: Fluid color with density-based alpha
                const [r, g, b] = params.fluidColor;
                const alpha = this.p.constrain(density * 3, 0, 255);

                this.buffer.pixels[pxIdx] = r;
                this.buffer.pixels[pxIdx + 1] = g;
                this.buffer.pixels[pxIdx + 2] = b;
                this.buffer.pixels[pxIdx + 3] = alpha;
            }
        }

        this.buffer.updatePixels();
        this.drawBuffer(params);
    }

    /**
     * Render fire simulation (temperature + smoke based)
     */
    renderFire(fire, params) {
        this.buffer.loadPixels();

        for (let y = 1; y <= this.N; y++) {
            for (let x = 1; x <= this.N; x++) {
                const intensity = fire.getIntensity(x, y);
                const pxIdx = 4 * ((x - 1) + (y - 1) * this.N);

                const [r, g, b, a] = this.getFireColor(intensity);

                this.buffer.pixels[pxIdx] = r;
                this.buffer.pixels[pxIdx + 1] = g;
                this.buffer.pixels[pxIdx + 2] = b;
                this.buffer.pixels[pxIdx + 3] = a;
            }
        }

        this.buffer.updatePixels();
        this.drawBuffer(params);
    }

    /**
     * Fire color mapping based on temperature and smoke
     */
    getFireColor(intensity) {
        const temp = intensity.temperature;
        const smoke = intensity.smoke;
        const reaction = intensity.reaction;

        // Normalize values with better scaling
        const t = Math.min(temp / 150, 1.0);      // Lower temp threshold
        const s = Math.min(smoke / 200, 1.0);      // Lower smoke threshold
        const r = Math.min(reaction * 5, 1.0);     // Reaction visibility

        let red, green, blue, alpha;

        // ACTIVE BURNING - Temperature-based gradient
        if (t > 0.8) {
            // White-hot core (very hot)
            red = 255;
            green = 240 + (t - 0.8) / 0.2 * 15;
            blue = 200 + (t - 0.8) / 0.2 * 55;
        } else if (t > 0.6) {
            // Yellow flame
            red = 255;
            green = 200 + (t - 0.6) / 0.2 * 40;
            blue = (t - 0.6) / 0.2 * 50;
        } else if (t > 0.4) {
            // Orange flame
            red = 255;
            green = 140 + (t - 0.4) / 0.2 * 60;
            blue = 0;
        } else if (t > 0.2) {
            // Red-orange
            red = 255;
            green = 50 + (t - 0.2) / 0.2 * 90;
            blue = 0;
        } else if (t > 0.08) {
            // Deep red / embers
            red = 180 + (t - 0.08) / 0.12 * 75;
            green = (t - 0.08) / 0.12 * 50;
            blue = 0;
        } else {
            // SMOKE (low temp, potentially high smoke)
            const smokeIntensity = s * 80;
            red = smokeIntensity;
            green = smokeIntensity;
            blue = smokeIntensity;
        }

        // Alpha: Visible if hot OR smoky
        // Increase visibility overall
        alpha = Math.min((t * 300) + (s * 150), 255);

        return [red, green, blue, alpha];
    }

    /**
     * Draw the buffer to canvas with blend mode
     */
    drawBuffer(params) {
        this.p.push();
        
        if (params.blendMode === 'additive') {
            this.p.blendMode(this.p.ADD);
        } else if (params.blendMode === 'multiply') {
            this.p.blendMode(this.p.MULTIPLY);
        }
        
        this.p.image(this.buffer, 0, 0, this.N * this.scale, this.N * this.scale);
        this.p.pop();
    }

    /**
     * Draw velocity vectors (works for both Fluid and Fire)
     */
    drawVectors(simulation, params) {
        const step = params.velocityVectorDensity;
        this.p.stroke(255, 150);
        this.p.strokeWeight(1);

        for (let y = 1; y <= this.N; y += step) {
            for (let x = 1; x <= this.N; x += step) {
                const i = simulation.idx(x, y);
                
                // Only draw where there's visible content
                let hasContent = false;
                if (simulation instanceof Fluid) {
                    hasContent = simulation.density[i] > 50;
                } else if (simulation instanceof Fire) {
                    hasContent = simulation.temperature[i] > 20 || simulation.smoke[i] > 50;
                }

                if (!hasContent) continue;

                const vx = simulation.vX[i] * params.velocityVectorScale * this.scale;
                const vy = simulation.vY[i] * params.velocityVectorScale * this.scale;

                if (vx * vx + vy * vy > 0.1) {
                    const sx = (x - 1) * this.scale;
                    const sy = (y - 1) * this.scale;
                    
                    this.p.line(sx, sy, sx + vx, sy + vy);
                    
                    // Arrowhead
                    this.p.push();
                    this.p.translate(sx + vx, sy + vy);
                    this.p.rotate(Math.atan2(vy, vx));
                    this.p.line(0, 0, -2, -2);
                    this.p.line(0, 0, -2, 2);
                    this.p.pop();
                }
            }
        }
    }
}