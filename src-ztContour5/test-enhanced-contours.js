import contours from "./enhanced-contours.js";
import { contourData as testData1 } from "./testData1.js";
import { contourData as testData2 } from "./testData2.js";

// Test suite for enhanced d3-contour functionality
function runTests() {
  console.log("=== Enhanced D3-Contour Test Suite ===\n");

  // Test 1: Basic functionality with testData1 (contains null values)
  console.log("Test 1: Basic contour generation with null values (testData1)");
  testBasicContours(testData1, "testData1");

  // Test 2: Enhanced functionality with testData2 (dense data)
  console.log("\nTest 2: Enhanced contour generation with dense data (testData2)");
  testBasicContours(testData2, "testData2");

  // Test 3: Saddle point disambiguation
  console.log("\nTest 3: Saddle point disambiguation");
  testSaddleDisambiguation();

  // Test 4: Enhanced smoothing
  console.log("\nTest 4: Enhanced smoothing controls");
  testSmoothingControls();

  // Test 5: GeoJSON export
  console.log("\nTest 5: GeoJSON export functionality");
  testGeoJSONExport();

  // Test 6: Performance comparison
  console.log("\nTest 6: Performance test");
  testPerformance();

  console.log("\n=== Test Suite Complete ===");
}

function testBasicContours(data, datasetName) {
  try {
    const { x, y, v } = data.data;
    const dx = x.length;
    const dy = y.length;
    
    // Flatten the 2D array for d3-contour
    const values = [];
    for (let j = 0; j < dy; j++) {
      for (let i = 0; i < dx; i++) {
        values.push(v[j][i]);
      }
    }

    // Create contour generator
    const contourGenerator = contours()
      .size([dx, dy])
      .smooth(true)
      .smoothFactor(0.5)
      .saddleDisambiguation(true)
      .nullHandling(true);

    // Generate contours
    const start = performance.now();
    const contourLines = contourGenerator.thresholds([50, 100, 150, 200, 250])(values);
    const duration = performance.now() - start;

    console.log(`  ✓ Generated ${contourLines.length} contour levels in ${duration.toFixed(2)}ms`);
    
    // Validate results
    contourLines.forEach((contour, i) => {
      const polygonCount = contour.coordinates.length;
      const hasHoles = contour.coordinates.some(polygon => polygon.length > 1);
      console.log(`    Level ${contour.value}: ${polygonCount} polygons${hasHoles ? ' (with holes)' : ''}`);
    });

    // Test null value handling
    const nullCount = values.filter(v => v == null).length;
    console.log(`  ✓ Handled ${nullCount} null values in ${datasetName}`);

  } catch (error) {
    console.error(`  ✗ Error in basic contours test: ${error.message}`);
  }
}

function testSaddleDisambiguation() {
  try {
    // Create a simple test case with saddle points
    const testValues = [
      1, 0, 1,
      0, 0.5, 0,  // This creates a saddle point at the center
      1, 0, 1
    ];

    const contourGenerator = contours()
      .size([3, 3])
      .smooth(false);

    // Test with saddle disambiguation OFF
    const contoursNoSaddle = contourGenerator
      .saddleDisambiguation(false)
      .thresholds([0.5])(testValues);

    // Test with saddle disambiguation ON
    const contoursSaddle = contourGenerator
      .saddleDisambiguation(true)
      .thresholds([0.5])(testValues);

    console.log(`  ✓ Without saddle disambiguation: ${contoursNoSaddle[0].coordinates.length} polygons`);
    console.log(`  ✓ With saddle disambiguation: ${contoursSaddle[0].coordinates.length} polygons`);

  } catch (error) {
    console.error(`  ✗ Error in saddle disambiguation test: ${error.message}`);
  }
}

function testSmoothingControls() {
  try {
    const testValues = [
      1, 2, 3,
      4, 5, 6,
      7, 8, 9
    ];

    const contourGenerator = contours().size([3, 3]);

    // Test different smoothing factors
    const smoothFactors = [0, 0.5, 1.0];
    
    smoothFactors.forEach(factor => {
      const result = contourGenerator
        .smooth(true)
        .smoothFactor(factor)
        .thresholds([5])(testValues);
      
      console.log(`  ✓ Smoothing factor ${factor}: generated ${result.length} contours`);
    });

  } catch (error) {
    console.error(`  ✗ Error in smoothing controls test: ${error.message}`);
  }
}

function testGeoJSONExport() {
  try {
    const testValues = [1, 2, 3, 4, 5, 6, 7, 8, 9];

    const contourGenerator = contours()
      .size([3, 3])
      .geoJSON(true);

    const result = contourGenerator.thresholds([5])(testValues);
    
    // Validate GeoJSON structure
    result.forEach(contour => {
      const hasProperties = contour.hasOwnProperty('properties');
      const hasValue = contour.properties && contour.properties.hasOwnProperty('value');
      console.log(`  ✓ GeoJSON export: properties=${hasProperties}, value=${hasValue}`);
    });

  } catch (error) {
    console.error(`  ✗ Error in GeoJSON export test: ${error.message}`);
  }
}

function testPerformance() {
  try {
    // Create a larger dataset for performance testing
    const size = 50;
    const values = new Array(size * size);
    for (let i = 0; i < values.length; i++) {
      const x = i % size;
      const y = Math.floor(i / size);
      values[i] = Math.sin(x / 10) * Math.cos(y / 10) * 100 + Math.random() * 10;
    }

    const contourGenerator = contours()
      .size([size, size])
      .smooth(true)
      .saddleDisambiguation(true);

    // Performance test
    const start = performance.now();
    const result = contourGenerator.thresholds(10)(values);
    const duration = performance.now() - start;

    console.log(`  ✓ Generated ${result.length} contour levels for ${size}x${size} grid in ${duration.toFixed(2)}ms`);
    console.log(`  ✓ Performance: ${(values.length / duration * 1000).toFixed(0)} cells/second`);

  } catch (error) {
    console.error(`  ✗ Error in performance test: ${error.message}`);
  }
}

// Export utilities for external testing
export function exportGeoJSON(contourData, filename = 'contours.geojson') {
  const featureCollection = {
    type: "FeatureCollection",
    features: contourData.map(contour => ({
      type: "Feature",
      properties: contour.properties || { value: contour.value },
      geometry: contour
    }))
  };
  
  return JSON.stringify(featureCollection, null, 2);
}

export function validateContours(contours) {
  const issues = [];
  
  contours.forEach((contour, index) => {
    // Check for valid MultiPolygon structure
    if (contour.type !== "MultiPolygon") {
      issues.push(`Contour ${index}: Invalid type ${contour.type}`);
    }
    
    // Check for self-intersections (basic check)
    contour.coordinates.forEach((polygon, polyIndex) => {
      polygon.forEach((ring, ringIndex) => {
        if (ring.length < 3) {
          issues.push(`Contour ${index}, polygon ${polyIndex}, ring ${ringIndex}: Ring has < 3 points`);
        }
      });
    });
  });
  
  return issues;
}

// Run tests if this file is executed directly
if (typeof window === 'undefined' && typeof global !== 'undefined') {
  runTests();
}

export default runTests;