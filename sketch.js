/**
 * sketch.js - The Controller
 * Manages UI, user input, and coordinates the Engine and Renderer.
 */

// ===== SIMULATION CONFIGURATION =====
const GRID_SIZE = 200; 
const CELL_SCALE = 3;  
const SOLVER_ITERATIONS = 20; 

let fluid;
let view;
let gui;

let prevMouseX = 0;
let prevMouseY = 0;

// ===== SIMULATION PARAMETERS (0-1 Sliders) =====
const parameters = {
    // Physics
    viscosity: 0.1,        
    diffusion: 0.0,
    dyeAmount: 0.4,        
    mouseForce: 0.2,       
    
    // Environmental
    buoyancy: 0.0,
    coolingRate: 0.1,      
    vorticityStrength: 0.2, 
    
    // Decay
    velocityDamping: 0.1,  
    densityFade: 0.05,     
    
    // Interaction
    brushRadius: 3,
    continuousInput: true,
    
    // Visual settings
    fluidColor: [200, 255, 255],
    backgroundColor: [10, 10, 15],
    blendMode: 'normal',
    showVelocityVectors: false,
    velocityVectorScale: 10,
    velocityVectorDensity: 8,
    
    resetSimulation: function() {
        const actualVisc = map(parameters.viscosity, 0, 1, 0, 0.000005);
        fluid = new Fluid(GRID_SIZE, 0.1, parameters.diffusion, actualVisc);
        console.log("Simulation reset");
    }
};

// Updated Presets (Smoke and Thick removed)
const presets = {
    fluid: { buoyancy: 0.0, coolingRate: 0.1, vorticityStrength: 0.1, densityFade: 0.05, velocityDamping: 0.1 },
    fire:  { buoyancy: 0.15, coolingRate: 0.4, vorticityStrength: 0.8, densityFade: 0.1, velocityDamping: 0.15 }
};

function setup() {
    createCanvas(GRID_SIZE * CELL_SCALE, GRID_SIZE * CELL_SCALE);
    
    const actualVisc = map(parameters.viscosity, 0, 1, 0, 0.000005);
    fluid = new Fluid(GRID_SIZE, 0.1, parameters.diffusion, actualVisc);
    view = new FluidRenderer(this, GRID_SIZE, CELL_SCALE);
    
    initializeGUI();
    prevMouseX = mouseX;
    prevMouseY = mouseY;
}

function draw() {
    background(parameters.backgroundColor[0], parameters.backgroundColor[1], parameters.backgroundColor[2]);

    if (mouseIsPressed) {
        processMouseInput();
    }

    // MATH MAPPING: 0-1 Sliders to Physics Values
    const physViscosity = map(parameters.viscosity, 0, 1, 0, 0.000005);
    const physVorticity = map(parameters.vorticityStrength, 0, 1, 0, 50);
    const physDamping = map(parameters.velocityDamping, 0, 1, 1.0, 0.9);
    const physFade = map(parameters.densityFade, 0, 1, 1.0, 0.9);
    const physCooling = map(parameters.coolingRate, 0, 1, 1.0, 0.9);

    fluid.viscosity = physViscosity;
    fluid.diffusion = parameters.diffusion;
    
    fluid.step({
        buoyancy: parameters.buoyancy,
        cooling: physCooling,
        damping: physDamping,
        fade: physFade,
        vorticity: physVorticity,
        iterations: SOLVER_ITERATIONS
    });

    view.render(fluid, parameters);
    
    if (parameters.showVelocityVectors) {
        view.drawVectors(fluid, parameters);
    }

    displayFrameRate();
    prevMouseX = mouseX;
    prevMouseY = mouseY;
}

function processMouseInput() {
    const gx = floor(mouseX / CELL_SCALE) + 1;
    const gy = floor(mouseY / CELL_SCALE) + 1;
    const mouseDeltaX = mouseX - prevMouseX;
    const mouseDeltaY = mouseY - prevMouseY;
    
    if (parameters.continuousInput || mouseDeltaX !== 0 || mouseDeltaY !== 0) {
        const r = parameters.brushRadius;
        const actualDye = map(parameters.dyeAmount, 0, 1, 50, 2000);
        const actualForce = map(parameters.mouseForce, 0, 1, 0.05, 2.0);

        for (let i = -r; i <= r; i++) {
            for (let j = -r; j <= r; j++) {
                if (i*i + j*j <= r*r) {
                    const weight = 1.0 - (Math.sqrt(i*i + j*j) / r);
                    const amt = actualDye * weight;
                    
                    fluid.addDensity(gx + i, gy + j, amt);
                    fluid.addVelocity(gx + i, gy + j, mouseDeltaX * actualForce, mouseDeltaY * actualForce);
                    
                    if (parameters.buoyancy > 0) {
                        fluid.addTemperature(gx + i, gy + j, amt * 0.8);
                    }
                }
            }
        }
    }
}

function initializeGUI() {
    gui = new dat.GUI();
    
    const presetsFolder = gui.addFolder('🎨 Presets');
    const presetController = { preset: 'fluid' };
    // This dynamically pulls only 'fluid' and 'fire' now
    presetsFolder.add(presetController, 'preset', Object.keys(presets)).name("Load Preset").onChange(applyPreset);
    presetsFolder.open();

    const physicsFolder = gui.addFolder('⚙️ Physics');
    physicsFolder.add(parameters, 'buoyancy', 0, 1).name("Buoyancy (Rise)");
    physicsFolder.add(parameters, 'vorticityStrength', 0, 1).name("Vorticity (Swirl)");
    physicsFolder.add(parameters, 'viscosity', 0, 1).name("Viscosity (Thick)");

    const decayFolder = gui.addFolder('⏱️ Decay');
    decayFolder.add(parameters, 'velocityDamping', 0, 1).name("Friction (Stop)").listen();
    decayFolder.add(parameters, 'densityFade', 0, 1).name("Fade (Ghost)").listen();
    decayFolder.add(parameters, 'coolingRate', 0, 1).name("Heat Cooling").listen();

    const inputFolder = gui.addFolder('🖱️ Interaction');
    inputFolder.add(parameters, 'brushRadius', 1, 15).name("Brush Size");
    inputFolder.add(parameters, 'dyeAmount', 0, 1).name("Dye Strength");
    inputFolder.add(parameters, 'mouseForce', 0, 1).name("Mouse Force");

    const visualFolder = gui.addFolder('🎨 Visuals');
    visualFolder.addColor(parameters, 'fluidColor').name("Fluid Color").listen();
    visualFolder.add(parameters, 'showVelocityVectors').name("Show Vectors");

    gui.add(parameters, 'resetSimulation').name("🔄 Reset Engine");
}

function applyPreset(name) {
    const preset = presets[name];
    for (let key in preset) {
        if (parameters.hasOwnProperty(key)) {
            parameters[key] = preset[key];
        }
    }
    for (let i in gui.__folders) {
        for (let j in gui.__folders[i].__controllers) {
            gui.__folders[i].__controllers[j].updateDisplay();
        }
    }
}

function displayFrameRate() {
    fill(255);
    noStroke();
    textSize(14);
    text(`FPS: ${floor(frameRate())} | Grid: ${GRID_SIZE}x${GRID_SIZE}`, 10, height - 10);
}