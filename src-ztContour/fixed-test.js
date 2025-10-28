/**
 * Test script for the fixed contour generator
 */

// Import the fixed contour implementation
import fixedContours, { toGeoJSON } from './fixed-contour.js';
import fs from 'fs';

// Create sample data
function createSampleData(width, height) {
    const data = new Array(width * height);
    
    // Create a test pattern with multiple peaks and valleys
    for (let y = 0; y < height; y++) {
        for (let x = 0; x < width; x++) {
            // Normalized coordinates
            const nx = x / width * 2 - 1;
            const ny = y / height * 2 - 1;
            
            // Multiple peaks function
            data[y * width + x] = 
                Math.exp(-8 * (nx * nx + ny * ny)) +  // Central peak
                Math.exp(-8 * ((nx - 0.5) * (nx - 0.5) + (ny - 0.5) * (ny - 0.5)) * 2) * 0.7 +  // Secondary peak
                Math.exp(-8 * ((nx + 0.5) * (nx + 0.5) + (ny + 0.5) * (ny + 0.5)) * 2) * 0.3;   // Tertiary peak
        }
    }
    
    return data;
}

// Generate contours for different modes
function runTest() {
    // Create data
    const width = 50;
    const height = 50;
    const data = createSampleData(width, height);
    
    // Create contour generator
    const contour = fixedContours()
        .size([width, height])
        .thresholds(10);
    
    // Generate surfaces
    contour.mode('surfaces');
    const surfaces = contour(data);
    fs.writeFileSync('testData-surfaces.json', JSON.stringify(toGeoJSON(surfaces), null, 2));
    
    // Generate lines
    contour.mode('lines');
    const lines = contour(data);
    fs.writeFileSync('testData-lines.json', JSON.stringify(toGeoJSON(lines), null, 2));
    
    // Generate both
    contour.mode('both');
    const both = contour(data);
    fs.writeFileSync('testData-both.json', JSON.stringify(toGeoJSON(both), null, 2));
    
    console.log('Generated test data files:');
    console.log('- testData-surfaces.json');
    console.log('- testData-lines.json');
    console.log('- testData-both.json');
}

// Run the test
runTest(); 