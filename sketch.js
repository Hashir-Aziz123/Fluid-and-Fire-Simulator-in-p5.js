/**
 * sketch.js - The Controller
 * Manages UI, mode switching, and coordinates simulations.
 */

// ===== SIMULATION CONFIGURATION =====
const GRID_SIZE = 200;
const CELL_SCALE = 3;
const SOLVER_ITERATIONS = 20;

let simulation; // Current active simulation (Fluid or Fire)
let renderer;
let gui;

let prevMouseX = 0;
let prevMouseY = 0;

// ===== MODE CONTROL =====
let currentMode = 'fluid'; // 'fluid' or 'fire'

// ===== SHARED PARAMETERS =====
const parameters = {
    // === MODE SELECTOR ===
    mode: 'fluid',
    
    // === PHYSICS (Shared) ===
    viscosity: 0.1,
    diffusion: 0.0,
    
    // === FLUID-SPECIFIC ===
    fluidDyeAmount: 0.4,
    fluidMouseForce: 0.2,
    fluidVorticity: 0.2,
    
    // === FIRE-SPECIFIC ===
    fireFuelAmount: 0.5,
    fireIgnitionTemp: 0.7,
    fireMouseForce: 0.15,
    fireVorticity: 0.6,
    
    // === DECAY (Different per mode) ===
    velocityDamping: 0.1,
    densityFade: 0.05,
    coolingRate: 0.1,
    buoyancy: 0.0,
    
    // === INTERACTION ===
    brushRadius: 3,
    continuousInput: true,
    
    // === VISUAL ===
    fluidColor: [200, 255, 255],
    backgroundColor: [10, 10, 15],
    blendMode: 'normal',
    showVelocityVectors: false,
    velocityVectorScale: 10,
    velocityVectorDensity: 8,
    
    // === ACTIONS ===
    resetSimulation: function() {
        createSimulation(parameters.mode);
        console.log(`${parameters.mode} simulation reset`);
    }
};

// ===== PRESETS =====
const presets = {
    fluid: {
        fluidDyeAmount: 0.4,
        fluidMouseForce: 0.2,
        fluidVorticity: 0.2,
        velocityDamping: 0.1,
        densityFade: 0.05,
        buoyancy: 0.0,
        coolingRate: 0.1
    },
   fire: {
        fireFuelAmount: 0.4,
        fireIgnitionTemp: 0.4,
        fireMouseForce: 0.1,
        fireVorticity: 0.6,    // Good swirl
        velocityDamping: 0.05, // Let it flow easier
        densityFade: 0.06,
        buoyancy: 0.25,        // Slightly less rise
        coolingRate: 0.2      // Higher cooling so it stops before top!
    }
};

function setup() {
    createCanvas(GRID_SIZE * CELL_SCALE, GRID_SIZE * CELL_SCALE);
    
    createSimulation(currentMode);
    renderer = new SimulationRenderer(this, GRID_SIZE, CELL_SCALE);
    
    initializeGUI();
    
    prevMouseX = mouseX;
    prevMouseY = mouseY;
}

function draw() {
    background(parameters.backgroundColor[0], 
               parameters.backgroundColor[1], 
               parameters.backgroundColor[2]);

    if (mouseIsPressed) {
        processMouseInput();
    }

    // Map UI sliders to physics values
    const physViscosity = map(parameters.viscosity, 0, 1, 0, 0.000005);
    const physDamping = map(parameters.velocityDamping, 0, 1, 1.0, 0.9);
    const physFade = map(parameters.densityFade, 0, 1, 1.0, 0.9);
    const physCooling = map(parameters.coolingRate, 0, 1, 1.0, 0.4);

    simulation.viscosity = physViscosity;
    simulation.diffusion = parameters.diffusion;

    // Step simulation
    if (simulation instanceof Fluid) {
        const physVorticity = map(parameters.fluidVorticity, 0, 1, 0, 50);
        
        simulation.step({
            vorticity: physVorticity,
            damping: physDamping,
            fade: physFade,
            iterations: SOLVER_ITERATIONS,
            buoyancy: parameters.buoyancy,
            cooling: physCooling
        });
        
        renderer.renderFluid(simulation, parameters);
    } else if (simulation instanceof Fire) {
       
        const fireSwirl = map(parameters.fireVorticity, 0, 1, 0, 10);

        simulation.step({
            buoyancy: parameters.buoyancy,
            damping: physDamping,
            fade: physFade,
            cooling: physCooling,
            iterations: SOLVER_ITERATIONS,
            vorticity: fireSwirl
        });
        
        renderer.renderFire(simulation, parameters);
    }

    if (parameters.showVelocityVectors) {
        renderer.drawVectors(simulation, parameters);
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

        if (simulation instanceof Fluid) {
            // FLUID MODE: Add dye and velocity
            const actualDye = map(parameters.fluidDyeAmount, 0, 1, 50, 2000);
            const actualForce = map(parameters.fluidMouseForce, 0, 1, 0.05, 2.0);

            for (let i = -r; i <= r; i++) {
                for (let j = -r; j <= r; j++) {
                    if (i*i + j*j <= r*r) {
                        const weight = 1.0 - (Math.sqrt(i*i + j*j) / r);
                        const amt = actualDye * weight;
                        
                        simulation.addDensity(gx + i, gy + j, amt);
                        simulation.addVelocity(gx + i, gy + j, 
                            mouseDeltaX * actualForce, 
                            mouseDeltaY * actualForce);
                    }
                }
            }
        } else if (simulation instanceof Fire) {
            // FIRE MODE: Add fuel, ignite, and velocity
            const actualFuel = map(parameters.fireFuelAmount, 0, 1, 100, 1000);
            const actualIgnition = map(parameters.fireIgnitionTemp, 0, 1, 30, 250);
            const actualForce = map(parameters.fireMouseForce, 0, 1, 0.05, 1.5);

            for (let i = -r; i <= r; i++) {
                for (let j = -r; j <= r; j++) {
                    if (i*i + j*j <= r*r) {
                        const dist = Math.sqrt(i*i + j*j);
                        const weight = 1.0 - (dist / r);
                        
                        // Add lots of fuel everywhere
                        simulation.addFuel(gx + i, gy + j, actualFuel * weight);
                        
                        // Add ignition heat broadly (not just center)
                        simulation.addTemperature(gx + i, gy + j, actualIgnition * weight * 0.7);
                        
                        // Add velocity
                        simulation.addVelocity(gx + i, gy + j, 
                            mouseDeltaX * actualForce, 
                            mouseDeltaY * actualForce);
                    }
                }
            }
        }
    }
}

function createSimulation(mode) {
    const actualVisc = map(parameters.viscosity, 0, 1, 0, 0.000005);
    
    if (mode === 'fluid') {
        simulation = new Fluid(GRID_SIZE, 0.1, parameters.diffusion, actualVisc);
    } else if (mode === 'fire') {
        simulation = new Fire(GRID_SIZE, 0.1, parameters.diffusion, actualVisc);
    }
    
    currentMode = mode;
}

function switchMode(newMode) {
    if (newMode !== currentMode) {
        createSimulation(newMode);
        applyPreset(newMode);
        
        // Adjust brush size based on mode
        if (newMode === 'fire') {
            parameters.brushRadius = 8;  // Larger brush for fire
        } else {
            parameters.brushRadius = 3;  // Standard brush for fluid
        }
        
        // Update GUI visibility
        updateGUIVisibility();
        
        // Update all GUI controllers to reflect new values
        for (let i in gui.__folders) {
            for (let j in gui.__folders[i].__controllers) {
                gui.__folders[i].__controllers[j].updateDisplay();
            }
        }
        
        console.log(`Switched to ${newMode} mode`);
    }
}

function initializeGUI() {
    gui = new dat.GUI();
    
    // === MODE SWITCHER ===
    const modeFolder = gui.addFolder('🔥💧 MODE');
    modeFolder.add(parameters, 'mode', ['fluid', 'fire'])
        .name("Simulation Type")
        .onChange(switchMode);
    modeFolder.open();

    // === PHYSICS ===
    const physicsFolder = gui.addFolder('⚙️ Physics');
    physicsFolder.add(parameters, 'viscosity', 0, 1).name("Viscosity");
    physicsFolder.add(parameters, 'buoyancy', 0, 1).name("Buoyancy (Rise)");
    
    // === FLUID CONTROLS ===
    window.fluidFolder = gui.addFolder('💧 Fluid Controls');
    fluidFolder.add(parameters, 'fluidVorticity', 0, 1).name("Vorticity (Swirl)");
    fluidFolder.add(parameters, 'fluidDyeAmount', 0, 1).name("Dye Amount");
    fluidFolder.add(parameters, 'fluidMouseForce', 0, 1).name("Mouse Force");
    
    // === FIRE CONTROLS ===
    window.fireFolder = gui.addFolder('🔥 Fire Controls');
    fireFolder.add(parameters, 'fireFuelAmount', 0, 1).name("Fuel Amount");
    fireFolder.add(parameters, 'fireIgnitionTemp', 0, 1).name("Ignition Heat");
    fireFolder.add(parameters, 'fireMouseForce', 0, 1).name("Mouse Force");
    fireFolder.add(parameters, 'fireVorticity', 0, 1).name("Vorticity (Turbulence)");

    // === DECAY ===
    const decayFolder = gui.addFolder('⏱️ Decay');
    decayFolder.add(parameters, 'velocityDamping', 0, 1).name("Velocity Damping");
    decayFolder.add(parameters, 'densityFade', 0, 1).name("Density Fade");
    decayFolder.add(parameters, 'coolingRate', 0, 1).name("Cooling Rate");

    // === INTERACTION ===
    const inputFolder = gui.addFolder('🖱️ Interaction');
    inputFolder.add(parameters, 'brushRadius', 1, 15).step(1).name("Brush Size");
    inputFolder.add(parameters, 'continuousInput').name("Continuous Input");

    // === VISUALS ===
    const visualFolder = gui.addFolder('🎨 Visuals');
    visualFolder.addColor(parameters, 'fluidColor').name("Fluid Color");
    visualFolder.addColor(parameters, 'backgroundColor').name("Background");
    visualFolder.add(parameters, 'showVelocityVectors').name("Show Vectors");

    // === RESET ===
    gui.add(parameters, 'resetSimulation').name("🔄 Reset Simulation");
    
    updateGUIVisibility();
}

function updateGUIVisibility() {
    if (currentMode === 'fluid') {
        fluidFolder.show();
        fireFolder.hide();
    } else if (currentMode === 'fire') {
        fluidFolder.hide();
        fireFolder.show();
    }
}

function applyPreset(name) {
    const preset = presets[name];
    if (!preset) return;
    
    for (let key in preset) {
        if (parameters.hasOwnProperty(key)) {
            parameters[key] = preset[key];
        }
    }
    
    // Update all GUI controllers
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
    text(`FPS: ${floor(frameRate())} | Mode: ${currentMode.toUpperCase()} | Grid: ${GRID_SIZE}x${GRID_SIZE}`, 
         10, height - 10);
}