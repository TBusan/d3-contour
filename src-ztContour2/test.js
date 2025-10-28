import fs from 'fs';
import { contours, contourDensity, toGeoJSON } from './index.js';
import { contourData } from './testData2.js';

// Function to generate test data
function generateTestData(width, height) {
  const values = new Float64Array(width * height);
  
  // Create a test pattern with multiple peaks and valleys
  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      // Normalized coordinates
      const nx = x / width * 2 - 1;
      const ny = y / height * 2 - 1;
      
      // Multiple peaks function
      values[y * width + x] = 
        Math.exp(-6 * (nx * nx + ny * ny)) +  // Central peak
        Math.exp(-8 * ((nx - 0.5) * (nx - 0.5) + (ny - 0.5) * (ny - 0.5)) * 2) * 0.7 +  // Secondary peak
        Math.exp(-8 * ((nx + 0.5) * (nx + 0.5) + (ny + 0.5) * (ny + 0.5)) * 2) * 0.5;   // Tertiary peak
    }
  }
  
  return values;
}

// Function to generate point data for density estimation
function generatePointData(numPoints) {
  const points = [];
  
  // Create clusters of points
  // Cluster 1: center
  for (let i = 0; i < numPoints * 0.4; i++) {
    points.push({
      x: Math.random() * 0.5 - 0.25 + 0.5,
      y: Math.random() * 0.5 - 0.25 + 0.5,
      value: Math.random() * 0.5 + 0.5
    });
  }
  
  // Cluster 2: top-right
  for (let i = 0; i < numPoints * 0.3; i++) {
    points.push({
      x: Math.random() * 0.3 + 0.7,
      y: Math.random() * 0.3 + 0.7,
      value: Math.random() * 0.5 + 0.5
    });
  }
  
  // Cluster 3: bottom-left
  for (let i = 0; i < numPoints * 0.3; i++) {
    points.push({
      x: Math.random() * 0.3 + 0.1,
      y: Math.random() * 0.3 + 0.1,
      value: Math.random() * 0.5 + 0.5
    });
  }
  
  return points;
}

// Test contour generation with grid data
function testContours() {
  const width = 50;
  const height = 50;
  const values = generateTestData(width, height);
  
  // Test contours with different modes
  console.log("Testing contour generation...");
  
  // Surfaces mode (default)
  const surfaceContours = contours()
    .size([width, height])
    .thresholds(10)
    .smooth(1.0)
    .mode("surfaces");
  
  const surfaces = surfaceContours(values);
  fs.writeFileSync('test-surfaces.json', JSON.stringify(toGeoJSON(surfaces), null, 2));
  
  // Lines mode
  const lineContours = contours()
    .size([width, height])
    .thresholds(10)
    .smooth(1.0)
    .mode("lines");
    
  const lines = lineContours(values);
  fs.writeFileSync('test-lines.json', JSON.stringify(toGeoJSON(lines), null, 2));
  
  // Both mode
  const bothContours = contours()
    .size([width, height])
    .thresholds(10)
    .smooth(1.0)
    .mode("both");
    
  const both = bothContours(values);
  fs.writeFileSync('test-both.json', JSON.stringify(toGeoJSON(both), null, 2));
  
  console.log("Contour tests completed. Output files:");
  console.log("- test-surfaces.json");
  console.log("- test-lines.json");
  console.log("- test-both.json");
}

// Test density estimation
function testDensity() {
  const points = generatePointData(1000);
  
  console.log("Testing density estimation...");
  
  // Create density estimator
  const density = contourDensity()
    .x(d => d.x)
    .y(d => d.y)
    .weight(d => d.value)
    .size([50, 50])
    .bandwidth(0.1)
    .thresholds(10)
    .kernelSize(5)
    .kernelDensityFactor(1.2)
    .mode("surfaces");
  
  const densityContours = density(points);
  fs.writeFileSync('test-density.json', JSON.stringify(toGeoJSON(densityContours), null, 2));
  
  console.log("Density test completed. Output file:");
  console.log("- test-density.json");
}

// Test contour generation with imported data
function testWithData1() {
  console.log("Testing contour generation with imported data...");
  
  const { data } = contourData;
  const { x, y, v } = data;
  
  // Use a much smaller subset of the data to avoid memory issues
  // Take every fourth point in both dimensions
  const subsampleFactor = 4;
  const subX = x.filter((_, i) => i % subsampleFactor === 0);
  const subY = y.filter((_, i) => i % subsampleFactor === 0);
  const subWidth = subX.length;
  const subHeight = subY.length;
  
  console.log(`Processing subsampled data: ${subWidth}x${subHeight} (original: ${x.length}x${y.length})`);
  
  // Convert 2D array to flat array for contour generation with subsampling
  const values = new Float64Array(subWidth * subHeight);
  for (let j = 0; j < subHeight; j++) {
    for (let i = 0; i < subWidth; i++) {
      const origJ = j * subsampleFactor;
      const origI = i * subsampleFactor;
      // Handle array bounds and null values
      if (origJ < v.length && origI < v[0].length) {
        values[j * subWidth + i] = v[origJ][origI];
      } else {
        values[j * subWidth + i] = null;
      }
    }
  }
  
  // Create contours with fewer thresholds
  const thresholdValues = [10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000];
  
  try {
    // Generate contours with both modes
    console.log("Generating contour lines...");
    const lineContours = contours()
      .size([subWidth, subHeight])
      .thresholds(thresholdValues)
      .smooth(1.0) // Disable smoothing to reduce computation
      .mode("lines")
      .nullValue(null); // Specify null value handling
      
    const lines = lineContours(values);
    fs.writeFileSync('test-data1-lines.json', JSON.stringify(toGeoJSON(lines), null, 2));
    
    console.log("Generating contour surfaces...");
    const surfaceContours = contours()
      .size([subWidth, subHeight])
      .thresholds(thresholdValues)
      .smooth(1.0) // Disable smoothing to reduce computation
      .mode("surfaces")
      .nullValue(null); // Specify null value handling
      
    const surfaces = surfaceContours(values);
    fs.writeFileSync('test-data1-surfaces.json', JSON.stringify(toGeoJSON(surfaces), null, 2));
    
    console.log("Data1 tests completed. Output files:");
    console.log("- test-data1-lines.json");
    console.log("- test-data1-surfaces.json");
  } catch (error) {
    console.error("Error processing contours:", error.message);
    console.error(error.stack);
  }
}

// Run tests
console.log("Running enhanced d3-contour tests...");
testContours();
testDensity();
testWithData1();
console.log("All tests completed successfully."); 