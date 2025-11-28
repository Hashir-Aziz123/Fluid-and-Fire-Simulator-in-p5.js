// --- CONFIGURATION ---
const N = 200;
const SCALE = 3;
const iter = 4;

let fluid;
let fluidImg;

// --- SETTINGS ---
let settings = {
    // Physics
    viscosity: 0.0000001,
    diffusion: 0.0000,
    dyeAmount: 250,
    force: 0.5,
    
    // Fire Physics
    buoyancy: 0.0,
    cooling: 0.99,
    turbulence: 0.5,
    
    // Interaction
    brushSize: 3,
    
    // Modes
    fireMode: false,
    
    // Visuals
    color: [200, 100, 255], // Changed to RGB format [R, G, B]
    showVelocity: false,
    vectorScale: 10,
    
    clear: function() {
        fluid = new Fluid(0.2, 0, 0.0000001);
    }
};

function setup() {
    let canvasSize = N * SCALE;
    createCanvas(canvasSize, canvasSize);
    noStroke();
    frameRate(60);
    
    // REMOVED HSB MODE - keep it in default RGB for dat.GUI compatibility
    // colorMode(HSB, 360, 100, 100);

    fluid = new Fluid(0.2, settings.diffusion, settings.viscosity);
    fluidImg = createImage(N, N); 

    // --- GUI ---
    let gui = new dat.GUI();
    
    // 1. Behavior Folder
    let f1 = gui.addFolder('Behavior');
    f1.add(settings, 'fireMode').name("🔥 Fire Mode").onChange(updateSimulationMode);
    
    f1.add(settings, 'buoyancy', 0, 0.10).name("Buoyancy").listen(); 
    f1.add(settings, 'cooling', 0.80, 1.0).name("Cooling").listen();
    f1.add(settings, 'turbulence', 0, 2.0).name("Wind/Chaos");
    f1.open();

    // 2. Interaction Folder
    let f2 = gui.addFolder('Interaction');
    f2.add(settings, 'dyeAmount', 100, 1000).name("Ink Amount");
    f2.add(settings, 'brushSize', 1, 10).name("Brush Size").step(1);
    f2.add(settings, 'force', 0.1, 3.0).name("Mouse Force");
    f2.open();

    // 3. Visuals Folder
    let f3 = gui.addFolder('Visuals');
    f3.addColor(settings, 'color').name("Fluid Color");
    f3.add(settings, 'showVelocity').name("Show Vectors");
    f3.add(settings, 'vectorScale', 0, 100).name("Arrow Size").listen(); 
    
    gui.add(settings, 'clear').name("Reset Fluid");
}

function updateSimulationMode() {
    if (settings.fireMode) {
        // FIRE PRESET
        settings.buoyancy = 0.08;
        settings.cooling = 0.90;
        settings.vectorScale = 10;
        settings.turbulence = 1.0;
    } else {
        // FLUID PRESET
        settings.buoyancy = 0.0;
        settings.cooling = 0.99;
        settings.vectorScale = 50;
        settings.turbulence = 0.0; // Fixed: was 0.5, should be 0.0 for fluid mode
    }
}

function draw() {
    background(0);

    // 1. SYNC SETTINGS TO ENGINE
    fluid.visc = settings.viscosity;
    fluid.diff = settings.diffusion;
    fluid.buoyancy = settings.buoyancy;
    fluid.cooling = settings.cooling;

    // 2. FIRE TURBULENCE
    if (settings.turbulence > 0) {
        addFireTurbulence();
    }

    // 3. INTERACTION & STEP
    handleMouse();
    fluid.step();

    // 4. RENDER FLUID
    renderFluid();

    // 5. RENDER VECTORS
    if (settings.showVelocity) {
        renderVelocityVectors();
    }
    
    fill(255);
    text("FPS: " + floor(frameRate()), 10, 20);
}

function handleMouse() {
    if (mouseIsPressed) {
        let centerX = floor(mouseX / SCALE);
        let centerY = floor(mouseY / SCALE);
        
        let radius = settings.brushSize;

        for (let i = -radius; i <= radius; i++) {
            for (let j = -radius; j <= radius; j++) {
                let x = centerX + i;
                let y = centerY + j;

                if (x > 0 && x < N - 1 && y > 0 && y < N - 1) {
                    if (i*i + j*j <= radius*radius) {
                        
                        let noiseAmt = random(0.5, 1.5);
                        
                        // A. Add Visuals (Dye)
                        fluid.addDensity(x, y, settings.dyeAmount * noiseAmt);
                        
                        // B. Add Physics (Heat)
                        if (settings.buoyancy > 0) {
                            fluid.addTemperature(x, y, 50 * noiseAmt);
                        }

                        // C. Add Motion (Velocity)
                        let amtX = movedX * settings.force;
                        let amtY = movedY * settings.force;
                        
                        if (settings.buoyancy > 0) {
                            let wiggle = random(-1, 1);
                            fluid.addVelocity(x, y, wiggle, -1);
                        }

                        fluid.addVelocity(x, y, amtX, amtY);
                    }
                }
            }
        }
    }
}

function addFireTurbulence() {
    let t = frameCount * 0.01;
    for (let j = 1; j < N - 1; j++) {
        for (let i = 1; i < N - 1; i++) {
            let index = IX(i, j);
            if (fluid.density[index] > 10) {
                let noiseVal = (noise(i * 0.1, j * 0.1, t) - 0.5); 
                fluid.Vx[index] += noiseVal * settings.turbulence;
            }
        }
    }
}

function renderFluid() {
    fluidImg.loadPixels();
    
    // FIXED: dat.GUI returns RGB array [r, g, b], use directly
    let rBase = settings.color[0];
    let gBase = settings.color[1];
    let bBase = settings.color[2];

    for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
            let d = fluid.density[IX(i, j)];
            let index = 4 * (i + j * N);
            
            if (d > 5) {
                if (settings.fireMode) {
                    // --- FIRE MODE ---
                    let r, g, b;
                    if (d > 200) { r=255; g=255; b=(d-200)*5; } 
                    else if (d > 100) { r=255; g=(d-100)*2.5; b=0; } 
                    else { r=d*2.5; g=d*0.5; b=0; } 
                    
                    fluidImg.pixels[index] = r;
                    fluidImg.pixels[index + 1] = g;
                    fluidImg.pixels[index + 2] = b;
                    fluidImg.pixels[index + 3] = 255;
                    
                } else {
                    // --- FLUID MODE (User Color) ---
                    // Use RGB values directly from dat.GUI
                    fluidImg.pixels[index]     = rBase;
                    fluidImg.pixels[index + 1] = gBase;
                    fluidImg.pixels[index + 2] = bBase;
                    
                    // Alpha mapping
                    fluidImg.pixels[index + 3] = constrain(d * 3, 0, 255); 
                }
            } else {
                fluidImg.pixels[index + 3] = 0;
            }
        }
    }
    fluidImg.updatePixels();
    image(fluidImg, 0, 0, width, height);
}

function renderVelocityVectors() {
    let step = 10;
    stroke(255, 150);
    strokeWeight(1);
    noFill();
    for (let i = 0; i < N; i += step) {
        for (let j = 0; j < N; j += step) {
            let x = i * SCALE;
            let y = j * SCALE;
            let index = IX(i, j);
            let vx = fluid.Vx[index];
            let vy = fluid.Vy[index];
            let len = Math.sqrt(vx*vx + vy*vy) * settings.vectorScale;
            
            if (len > 1) {
                len = constrain(len, 2, 40);
                let angle = atan2(vy, vx);
                push();
                translate(x, y);
                rotate(angle);
                line(0, 0, len, 0);
                let arrowSize = 3;
                line(len, 0, len - arrowSize, -arrowSize);
                line(len, 0, len - arrowSize, arrowSize);
                pop();
            }
        }
    }
}

function IX(x, y) {
    x = constrain(x, 0, N + 1);
    y = constrain(y, 0, N + 1);
    return x + (N + 2) * y;
}