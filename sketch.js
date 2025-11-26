// --- CONFIGURATION ---
const N = 128;         // Grid Resolution
const SCALE = 5;       // Scale factor
const iter = 4;        // Solver iterations (used by fluid.js)

// The Fluid Object
let fluid;

function setup() {
    let canvasSize = N * SCALE;
    createCanvas(canvasSize, canvasSize);
    noStroke();
    frameRate(60);

    // Create the fluid simulation
    // Fluid(dt, diffusion, viscosity)
    fluid = new Fluid(0.2, 0, 0.0000001);
}

function draw() {
    background(0);

    // 1. Interactions
    handleMouse();

    // 2. Physics Update
    fluid.step();

    // 3. Render
    // We access fluid.density directly to draw it
    for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
            let x = i * SCALE;
            let y = j * SCALE;
            
            // Get density from fluid object
            let d = fluid.density[IX(i, j)];

            if (d > 0.1) {
                fill(d);
                rect(x, y, SCALE, SCALE);
            }
        }
    }
    
    // Debug info
    fill(255);
    text("FPS: " + floor(frameRate()), 10, 20);
}

function handleMouse() {
    if (mouseIsPressed) {
        let x = floor(mouseX / SCALE);
        let y = floor(mouseY / SCALE);
        
        // Add Dye
        fluid.addDensity(x, y, 200);
        
        // Add Velocity
        let amtX = movedX * 0.5;
        let amtY = movedY * 0.5;
        fluid.addVelocity(x, y, amtX, amtY);
    }
}