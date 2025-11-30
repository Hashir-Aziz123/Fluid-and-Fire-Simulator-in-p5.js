// ===== SIMULATION CONFIGURATION =====
const GRID_SIZE = 200;
const CELL_SCALE = 3;
const SOLVER_ITERATIONS = 4;

let fluidSimulation;
let renderBuffer;
let gui;

// ===== SIMULATION PARAMETERS =====
const parameters = {
    // Fluid physics properties (CRITICAL: These must stay very small!)
    viscosity: 0.0000001,
    diffusion: 0.0,
    dyeAmount: 250,
    mouseForce: 0.5,
    
    // Fire-specific physics
    buoyancy: 0.0,
    coolingRate: 0.99,
    turbulenceStrength: 0.0,
    
    // Decay rates (Keep these high for lingering fluid!)
    velocityDamping: 0.99,
    densityFade: 0.995,
    
    // User interaction
    brushRadius: 3,
    continuousInput: true,
    
    // Simulation modes
    isFireMode: false,
    
    // Visual settings
    fluidColor: [200, 255, 255],
    backgroundColor: [0, 0, 0],
    blendMode: 'normal',
    fadeTrails: false,
    showVelocityVectors: false,
    velocityVectorScale: 10,
    velocityVectorDensity: 10,
    
    // Performance
    simulationSpeed: 1.0,
    
    // Presets
    loadPreset: function(presetName) {
        applyPreset(presetName);
    },
    
    // Utility functions
    resetSimulation: function() {
        fluidSimulation = new Fluid(0.2, parameters.diffusion, parameters.viscosity);
    }
};

// Preset configurations (FIXED: Actual fluid behavior, not brush tool)
const presets = {
    fluid: {
        viscosity: 0.0000001,
        diffusion: 0.0,
        buoyancy: 0.0,
        coolingRate: 0.99,
        turbulenceStrength: 0.0,
        velocityDamping: 0.99,
        densityFade: 0.995,
        dyeAmount: 250,
        mouseForce: 0.5,
        fluidColor: [200, 255, 255],
        isFireMode: false
    },
    fire: {
        viscosity: 0.0000001,
        diffusion: 0.0,
        buoyancy: 0.08,
        coolingRate: 0.90,
        turbulenceStrength: 1.0,
        velocityDamping: 0.99,
        densityFade: 0.995,
        dyeAmount: 250,
        mouseForce: 0.5,
        isFireMode: true
    },
    smoke: {
        viscosity: 0.0000001,
        diffusion: 0.0,
        buoyancy: 0.03,
        coolingRate: 0.97,
        turbulenceStrength: 0.3,
        velocityDamping: 0.99,
        densityFade: 0.992,
        dyeAmount: 300,
        mouseForce: 0.5,
        fluidColor: [180, 180, 180],
        isFireMode: false
    },
    thick: {
        viscosity: 0.000001,
        diffusion: 0.0,
        buoyancy: 0.0,
        coolingRate: 0.99,
        turbulenceStrength: 0.0,
        velocityDamping: 0.97,
        densityFade: 0.998,
        dyeAmount: 400,
        mouseForce: 0.3,
        fluidColor: [100, 200, 255],
        isFireMode: false
    },
    wispy: {
        viscosity: 0.00000001,
        diffusion: 0.00001,
        buoyancy: 0.01,
        coolingRate: 0.98,
        turbulenceStrength: 0.5,
        velocityDamping: 0.985,
        densityFade: 0.99,
        dyeAmount: 200,
        mouseForce: 0.7,
        fluidColor: [255, 200, 255],
        isFireMode: false
    }
};

function setup() {
    const canvasSize = GRID_SIZE * CELL_SCALE;
    createCanvas(canvasSize, canvasSize);
    noStroke();
    frameRate(60);
    
    fluidSimulation = new Fluid(0.2, parameters.diffusion, parameters.viscosity);
    renderBuffer = createImage(GRID_SIZE, GRID_SIZE);
    
    initializeGUI();
}

function initializeGUI() {
    gui = new dat.GUI();
    
    // Presets folder
    const presetsFolder = gui.addFolder('🎨 Presets');
    const presetController = { preset: 'fluid' };
    presetsFolder.add(presetController, 'preset', ['fluid', 'fire', 'smoke', 'thick', 'wispy'])
        .name("Load Preset")
        .onChange((value) => {
            applyPreset(value);
        });
    presetsFolder.open();
    
    // Behavior controls
    const behaviorFolder = gui.addFolder('⚙️ Physics (Keep values SMALL!)');
    behaviorFolder.add(parameters, 'isFireMode')
        .name("Fire Mode")
        .onChange(switchSimulationMode);
    behaviorFolder.add(parameters, 'buoyancy', 0, 0.15)
        .name("Buoyancy (upward force)")
        .step(0.001)
        .listen();
    behaviorFolder.add(parameters, 'coolingRate', 0.80, 1.0)
        .name("Cooling (0.9=fast, 0.99=slow)")
        .step(0.001)
        .listen();
    behaviorFolder.add(parameters, 'turbulenceStrength', 0, 2.0)
        .name("Turbulence")
        .step(0.1)
        .listen();
    behaviorFolder.add(parameters, 'viscosity', 0, 0.000005)
        .name("Viscosity (thickness)")
        .step(0.0000001)
        .listen();
    behaviorFolder.add(parameters, 'diffusion', 0, 0.0001)
        .name("Diffusion (spreading)")
        .step(0.000001)
        .listen();

    // Decay controls
    const decayFolder = gui.addFolder('⏱️ Decay (higher = lingers longer)');
    decayFolder.add(parameters, 'velocityDamping', 0.95, 0.999)
        .name("Velocity Damping")
        .step(0.001);
    decayFolder.add(parameters, 'densityFade', 0.98, 0.999)
        .name("Density Fade")
        .step(0.001);

    // Interaction controls
    const interactionFolder = gui.addFolder('🖱️ Interaction');
    interactionFolder.add(parameters, 'dyeAmount', 50, 1000)
        .name("Dye Amount")
        .step(10)
        .listen();
    interactionFolder.add(parameters, 'brushRadius', 1, 15)
        .name("Brush Size")
        .step(1);
    interactionFolder.add(parameters, 'mouseForce', 0.1, 5.0)
        .name("Mouse Force")
        .step(0.1)
        .listen();
    interactionFolder.add(parameters, 'continuousInput')
        .name("Continuous Input");
    interactionFolder.open();

    // Visual controls
    const visualFolder = gui.addFolder('🎨 Visuals');
    visualFolder.addColor(parameters, 'fluidColor')
        .name("Fluid Color")
        .listen();
    visualFolder.addColor(parameters, 'backgroundColor')
        .name("Background");
    visualFolder.add(parameters, 'blendMode', ['normal', 'additive', 'multiply'])
        .name("Blend Mode");
    visualFolder.add(parameters, 'fadeTrails')
        .name("Fade Trails (motion blur)");
    visualFolder.add(parameters, 'showVelocityVectors')
        .name("Show Vectors");
    visualFolder.add(parameters, 'velocityVectorScale', 0, 100)
        .name("Vector Scale")
        .step(1)
        .listen();
    visualFolder.add(parameters, 'velocityVectorDensity', 5, 20)
        .name("Vector Density")
        .step(1);
    
    // Performance controls
    const performanceFolder = gui.addFolder('⚡ Performance');
    performanceFolder.add(parameters, 'simulationSpeed', 0.1, 2.0)
        .name("Simulation Speed")
        .step(0.1);
    
    gui.add(parameters, 'resetSimulation').name("🔄 Reset");
}

function applyPreset(presetName) {
    const preset = presets[presetName];
    if (!preset) return;
    
    // Apply all preset values
    Object.keys(preset).forEach(key => {
        if (parameters.hasOwnProperty(key)) {
            parameters[key] = preset[key];
        }
    });
    
    // Update GUI displays
    for (let folderName in gui.__folders) {
        gui.__folders[folderName].__controllers.forEach(controller => {
            controller.updateDisplay();
        });
    }
    
    // Sync to simulation
    syncParametersToSimulation();
}

function switchSimulationMode() {
    if (parameters.isFireMode) {
        // Fire mode: hot, fast-moving, chaotic
        parameters.buoyancy = 0.08;
        parameters.coolingRate = 0.90;
        parameters.velocityVectorScale = 10;
        parameters.turbulenceStrength = 1.0;
    } else {
        // Fluid mode: calm, lingering, smooth
        parameters.buoyancy = 0.0;
        parameters.coolingRate = 0.99;
        parameters.velocityVectorScale = 50;
        parameters.turbulenceStrength = 0.0;
    }
    
    // Force GUI to update all controllers in all folders
    for (let folderName in gui.__folders) {
        gui.__folders[folderName].__controllers.forEach(controller => {
            controller.updateDisplay();
        });
    }
}

function draw() {
    // Apply background with optional fade trails
    if (parameters.fadeTrails) {
        fill(parameters.backgroundColor[0], parameters.backgroundColor[1], parameters.backgroundColor[2], 10);
        rect(0, 0, width, height);
    } else {
        background(parameters.backgroundColor[0], parameters.backgroundColor[1], parameters.backgroundColor[2]);
    }
    
    syncParametersToSimulation();
    
    if (parameters.turbulenceStrength > 0) {
        applyTurbulence();
    }
    
    processMouseInput();
    
    // Apply simulation speed multiplier
    for (let i = 0; i < parameters.simulationSpeed; i++) {
        fluidSimulation.step();
    }
    
    renderFluidDensity();
    
    if (parameters.showVelocityVectors) {
        renderVelocityField();
    }
    
    displayFrameRate();
}

function syncParametersToSimulation() {
    fluidSimulation.viscosity = parameters.viscosity;
    fluidSimulation.diffusion = parameters.diffusion;
    fluidSimulation.buoyancy = parameters.buoyancy;
    fluidSimulation.coolingRate = parameters.coolingRate;
    fluidSimulation.velocityDamping = parameters.velocityDamping;
    fluidSimulation.densityFade = parameters.densityFade;
}

function processMouseInput() {
    const shouldProcess = parameters.continuousInput ? mouseIsPressed : mouseIsPressed && (movedX !== 0 || movedY !== 0);
    
    if (!shouldProcess) return;
    
    const gridX = floor(mouseX / CELL_SCALE);
    const gridY = floor(mouseY / CELL_SCALE);
    const radius = parameters.brushRadius;
    
    // Apply forces in a circular brush pattern
    for (let offsetX = -radius; offsetX <= radius; offsetX++) {
        for (let offsetY = -radius; offsetY <= radius; offsetY++) {
            const cellX = gridX + offsetX;
            const cellY = gridY + offsetY;
            
            // Check grid bounds and circular shape
            const isInsideGrid = cellX > 0 && cellX < GRID_SIZE - 1 && 
                                cellY > 0 && cellY < GRID_SIZE - 1;
            const distanceSquared = offsetX * offsetX + offsetY * offsetY;
            const isInsideCircle = distanceSquared <= radius * radius;
            
            if (isInsideGrid && isInsideCircle) {
                // Softer falloff from center
                const falloff = 1.0 - (Math.sqrt(distanceSquared) / radius);
                const intensityVariation = random(0.7, 1.3) * falloff;
                
                // Add visual dye/smoke
                fluidSimulation.addDensity(
                    cellX, 
                    cellY, 
                    parameters.dyeAmount * intensityVariation
                );
                
                // Add heat (only when buoyancy is active)
                if (parameters.buoyancy > 0) {
                    fluidSimulation.addTemperature(
                        cellX, 
                        cellY, 
                        50 * intensityVariation
                    );
                }
                
                // Add velocity based on mouse movement
                const velocityX = movedX * parameters.mouseForce;
                const velocityY = movedY * parameters.mouseForce;
                
                // Add upward push for fire/buoyancy
                if (parameters.buoyancy > 0) {
                    const horizontalWiggle = random(-0.5, 0.5);
                    fluidSimulation.addVelocity(cellX, cellY, horizontalWiggle, -0.8);
                }
                
                fluidSimulation.addVelocity(cellX, cellY, velocityX, velocityY);
            }
        }
    }
}

function applyTurbulence() {
    const timeOffset = frameCount * 0.01;
    
    for (let y = 1; y < GRID_SIZE - 1; y++) {
        for (let x = 1; x < GRID_SIZE - 1; x++) {
            const index = getGridIndex(x, y);
            
            // Only apply wind to areas with existing density
            if (fluidSimulation.density[index] > 10) {
                const noiseValue = (noise(x * 0.1, y * 0.1, timeOffset) - 0.5);
                fluidSimulation.horizontalVelocity[index] += 
                    noiseValue * parameters.turbulenceStrength;
            }
        }
    }
}

function renderFluidDensity() {
    renderBuffer.loadPixels();
    
    const redChannel = parameters.fluidColor[0];
    const greenChannel = parameters.fluidColor[1];
    const blueChannel = parameters.fluidColor[2];
    
    for (let x = 0; x < GRID_SIZE; x++) {
        for (let y = 0; y < GRID_SIZE; y++) {
            const densityValue = fluidSimulation.density[getGridIndex(x, y)];
            const pixelIndex = 4 * (x + y * GRID_SIZE);
            
            if (densityValue > 5) {
                if (parameters.isFireMode) {
                    // Fire has dynamic color based on temperature
                    let red, green, blue;
                    
                    if (densityValue > 200) {
                        // White-hot core
                        red = 255;
                        green = 255;
                        blue = constrain((densityValue - 200) * 5, 0, 255);
                    } else if (densityValue > 100) {
                        // Yellow-orange flames
                        red = 255;
                        green = constrain((densityValue - 100) * 2.5, 0, 255);
                        blue = 0;
                    } else {
                        // Red embers
                        red = constrain(densityValue * 2.5, 0, 255);
                        green = constrain(densityValue * 0.5, 0, 255);
                        blue = 0;
                    }
                    
                    renderBuffer.pixels[pixelIndex] = red;
                    renderBuffer.pixels[pixelIndex + 1] = green;
                    renderBuffer.pixels[pixelIndex + 2] = blue;
                    renderBuffer.pixels[pixelIndex + 3] = 255;
                } else {
                    // Fluid mode uses user-selected color with density-based opacity
                    renderBuffer.pixels[pixelIndex] = redChannel;
                    renderBuffer.pixels[pixelIndex + 1] = greenChannel;
                    renderBuffer.pixels[pixelIndex + 2] = blueChannel;
                    renderBuffer.pixels[pixelIndex + 3] = constrain(densityValue * 3, 0, 255);
                }
            } else {
                // Transparent for low density areas
                renderBuffer.pixels[pixelIndex + 3] = 0;
            }
        }
    }
    
    renderBuffer.updatePixels();
    
    // Apply blend mode
    push();
    if (parameters.blendMode === 'additive') {
        blendMode(ADD);
    } else if (parameters.blendMode === 'multiply') {
        blendMode(MULTIPLY);
    }
    image(renderBuffer, 0, 0, width, height);
    pop();
}

function renderVelocityField() {
    const samplingStep = parameters.velocityVectorDensity;
    
    stroke(255, 150);
    strokeWeight(1);
    noFill();
    
    for (let x = 0; x < GRID_SIZE; x += samplingStep) {
        for (let y = 0; y < GRID_SIZE; y += samplingStep) {
            const screenX = x * CELL_SCALE;
            const screenY = y * CELL_SCALE;
            const index = getGridIndex(x, y);
            
            const velocityX = fluidSimulation.horizontalVelocity[index];
            const velocityY = fluidSimulation.verticalVelocity[index];
            const magnitude = Math.sqrt(velocityX * velocityX + velocityY * velocityY);
            const scaledLength = magnitude * parameters.velocityVectorScale;
            
            if (scaledLength > 1) {
                const arrowLength = constrain(scaledLength, 2, 40);
                const angle = atan2(velocityY, velocityX);
                const arrowheadSize = 3;
                
                push();
                translate(screenX, screenY);
                rotate(angle);
                
                // Draw arrow shaft
                line(0, 0, arrowLength, 0);
                
                // Draw arrowhead
                line(arrowLength, 0, arrowLength - arrowheadSize, -arrowheadSize);
                line(arrowLength, 0, arrowLength - arrowheadSize, arrowheadSize);
                
                pop();
            }
        }
    }
}

function displayFrameRate() {
    fill(255);
    noStroke();
    text("FPS: " + floor(frameRate()), 10, 20);
}

function getGridIndex(x, y) {
    x = constrain(x, 0, GRID_SIZE + 1);
    y = constrain(y, 0, GRID_SIZE + 1);
    return x + (GRID_SIZE + 2) * y;
}