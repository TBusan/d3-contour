/**
 * Test file for the enhanced contour library
 * Shows how to generate contours from both grid data and point data
 */

import { contours, contourDensity, toGeoJSON } from './index.js';

// Function to generate sample data
function generateSampleData(size, generator) {
    const [width, height] = size;
    const data = new Float64Array(width * height);
    
    for (let y = 0; y < height; y++) {
        for (let x = 0; x < width; x++) {
            data[y * width + x] = generator(x, y, width, height);
        }
    }
    
    return data;
}

// Generate a simple peak function
function samplePeak(x, y, width, height) {
    const centerX = width / 2;
    const centerY = height / 2;
    const dx = x - centerX;
    const dy = y - centerY;
    const distance = Math.sqrt(dx * dx + dy * dy);
    const peak = 100 * Math.exp(-distance * distance / (width * height / 10));
    
    // Add some noise
    return peak + Math.random() * 5;
}

// Test basic contours on grid data
function testBasicContours() {
    console.log("Testing basic contours...");
    
    // Generate sample data
    const size = [100, 100];
    const data = generateSampleData(size, samplePeak);
    
    // Create contour generator with optimized settings
    const contour = contours()
        .size(size)
        .thresholds(10)         // Generate 10 contour levels
        .smooth(0.5)            // Reduced smoothing level
        .mode('both');          // Generate both lines and surfaces
    
    // Generate contours
    const result = contour(data);
    
    // Convert to GeoJSON for use with mapping libraries
    const geoJson = toGeoJSON(result);
    
    console.log(`Generated ${result.length} contour levels`);
    console.log("First contour level value:", result[0].value);
    console.log("GeoJSON type:", geoJson.type);
    
    return result;
}

// Test density contours on point data
function testDensityContours() {
    console.log("Testing density contours...");
    
    // Generate some point data
    const pointCount = 1000;
    const points = [];
    const width = 500;
    const height = 500;
    
    for (let i = 0; i < pointCount; i++) {
        const x = width / 2 + (Math.random() - 0.5) * width * 0.8;
        const y = height / 2 + (Math.random() - 0.5) * height * 0.8;
        points.push({
            x: x,
            y: y,
            value: Math.random() * 10
        });
    }
    
    // Create density contour generator
    const density = contourDensity()
        .x(d => d.x)
        .y(d => d.y)
        .weight(d => d.value)
        .size([width, height])
        .bandwidth(20)          // Adjust bandwidth for smoothing
        .thresholds(10)         // Number of contour levels
        .bandwidthAdjust(1.5)   // Increase bandwidth for smoother contours
        .mode('both');          // Generate both lines and surfaces
    
    // Generate density contours
    const result = density(points);
    
    console.log(`Generated ${result.length} density contour levels`);
    if (result.length > 0) {
        console.log("First density contour level value:", result[0].value);
    } else {
        console.log("No density contours generated");
    }
    
    return result;
}

// Run tests
const basicResult = testBasicContours();
const densityResult = testDensityContours();

export { basicResult, densityResult };
