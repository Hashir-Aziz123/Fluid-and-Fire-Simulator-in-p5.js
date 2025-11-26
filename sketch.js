// --- CONFIGURATION ---
const N = 128;
const SCALE = 5;
const iter = 4;

let fluid;

// --- GUI SETTINGS ---
let settings = {
    viscosity: 0.0000001,
    diffusion: 0.0000,
    dyeAmount: 200, 
    force: 0.5,
    
    // Visuals
    fireMode: true,       // Default to Fire Mode
    color: [0, 100, 100], // Fallback color (if fire mode is off)
    
    clear: function() {
        fluid = new Fluid(0.2, 0, 0.0000001);
    }
};

function setup() {
    let canvasSize = N * SCALE;
    createCanvas(canvasSize, canvasSize);
    noStroke();
    frameRate(60);
    
    // HSB Mode: Hue(0-360), Saturation(0-100), Brightness(0-100), Alpha(0-1)
    colorMode(HSB, 360, 100, 100, 1);

    fluid = new Fluid(0.2, settings.diffusion, settings.viscosity);

    // --- SETUP GUI ---
    let gui = new dat.GUI();
    gui.add(settings, 'viscosity', 0, 0.005).name("Viscosity").step(0.000001);
    gui.add(settings, 'diffusion', 0, 0.01).name("Diffusion").step(0.000001);
    gui.add(settings, 'force', 0.1, 2.0).name("Mouse Force");
    gui.add(settings, 'dyeAmount', 100, 2000).name("Dye Amount");
    
    let f1 = gui.addFolder('Visuals');
    f1.add(settings, 'fireMode').name("🔥 Fire Mode"); // New Checkbox
    f1.addColor(settings, 'color').name("Custom Color");
    f1.open();
    
    gui.add(settings, 'clear').name("Reset Fluid");
}

function draw() {
    background(0); // Black background

    // Update settings
    fluid.visc = settings.viscosity;
    fluid.diff = settings.diffusion;

    handleMouse();
    fluid.step();

    // --- RENDER LOOP ---
    for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
            let x = i * SCALE;
            let y = j * SCALE;
            
            let d = fluid.density[IX(i, j)];

            if (d > 2) { 
                if (settings.fireMode) {
                    // --- FIRE MODE LOGIC ---
                    // Map density 'd' to specific fire colors
                    // d goes from 0 to ~255 (or higher)
                    
                    let h, s, b;
                    
                    // 1. Brightness is simply the density (capped at 100)
                    b = constrain(d, 0, 100);
                    
                    if (d > 200) {
                        // CORE: White (Too hot for color)
                        h = 0; 
                        s = 0; // Desaturate to white
                    } else if (d > 120) {
                        // INNER: Blue (Gas flame core)
                        h = 220; 
                        s = 100;
                    } else if (d > 60) {
                        // MID: Orange/Yellow
                        h = 30; 
                        s = 100;
                    } else {
                        // OUTER: Red (Cooling down)
                        h = 0; 
                        s = 100;
                    }
                    
                    fill(h, s, b);
                    
                } else {
                    // --- CUSTOM SINGLE COLOR MODE ---
                    let h = settings.color[0];
                    let s = settings.color[1];
                    // Fade out brightness based on density
                    let b = constrain(d, 0, settings.color[2]); 
                    fill(h, s, b);
                }
                
                rect(x, y, SCALE, SCALE);
            }
        }
    }
    
    fill(255);
    text("FPS: " + floor(frameRate()), 10, 20);
}

function handleMouse() {
    if (mouseIsPressed) {
        let x = floor(mouseX / SCALE);
        let y = floor(mouseY / SCALE);
        
        fluid.addDensity(x, y, settings.dyeAmount);
        
        let amtX = movedX * settings.force;
        let amtY = movedY * settings.force;
        fluid.addVelocity(x, y, amtX, amtY);
    }
}

function IX(x, y) {
    x = constrain(x, 0, N + 1);
    y = constrain(y, 0, N + 1);
    return x + (N + 2) * y;
}