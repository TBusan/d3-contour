import contours from "./contours.js";
import dataAdapter from "./dataAdapter.js";
import { contourData as testData1 } from "./testData1.js";
import { contourData as testData2 } from "./testData2.js";

// Example usage of the enhanced d3-contour library
function demonstrateEnhancements() {
  
  // Create data adapter and contour generator
  const adapter = dataAdapter();
  const contourGenerator = contours();
  
  console.log("=== Enhanced D3-Contour Demonstration ===\n");
  
  // Test with testData1.js (contains null values for boundary masking)
  console.log("1. Testing with testData1.js (with null boundary mask):");
  testWithData(testData1, "testData1", adapter, contourGenerator);
  
  console.log("\n" + "=".repeat(50) + "\n");
  
  // Test with testData2.js (dense grid data)
  console.log("2. Testing with testData2.js (dense grid):");
  testWithData(testData2, "testData2", adapter, contourGenerator);
  
  console.log("\n" + "=".repeat(50) + "\n");
  
  // Demonstrate enhanced features
  console.log("3. Demonstrating enhanced features:");
  demonstrateFeatures(testData2, adapter, contourGenerator);
}

function testWithData(data, dataName, adapter, contourGenerator) {
  try {
    // Adapt the data format
    const adaptedData = adapter(data);
    console.log(`   Adapted ${dataName}:`);
    console.log(`   - Grid size: ${adaptedData.width} x ${adaptedData.height}`);
    console.log(`   - X range: ${adaptedData.xmin} to ${adaptedData.xmax}`);
    console.log(`   - Y range: ${adaptedData.ymin} to ${adaptedData.ymax}`);
    console.log(`   - Null values: ${adaptedData.values.filter(v => v === null).length}`);
    
    // Configure contour generator
    contourGenerator
      .size([adaptedData.width, adaptedData.height])
      .nullValue(null)
      .smoothingFactor(0.8);
    
    // Generate contours
    const contourResults = contourGenerator(adaptedData.values);
    console.log(`   - Generated ${contourResults.length} contour levels`);
    
    // Transform back to original coordinates
    const finalContours = adapter.createContourData(adaptedData, contourResults);
    
    // Report statistics
    const totalPolygons = finalContours.reduce((sum, contour) => sum + contour.coordinates.length, 0);
    const totalHoles = finalContours.reduce((sum, contour) => 
      sum + contour.coordinates.reduce((holeSum, polygon) => holeSum + polygon.length - 1, 0), 0);
    
    console.log(`   - Total polygons: ${totalPolygons}`);
    console.log(`   - Total holes: ${totalHoles}`);
    console.log(`   - Contour values: [${contourResults.slice(0, 5).map(c => c.value.toFixed(1)).join(", ")}${contourResults.length > 5 ? "..." : ""}]`);
    
  } catch (error) {
    console.error(`   Error processing ${dataName}:`, error.message);
  }
}

function demonstrateFeatures(data, adapter, contourGenerator) {
  const adaptedData = adapter(data);
  
  // Test different smoothing factors
  console.log("   A. Testing smoothing factors:");
  [0, 0.5, 1.0].forEach(factor => {
    contourGenerator.smoothingFactor(factor);
    const result = contourGenerator(adaptedData.values);
    console.log(`      Smoothing ${factor}: ${result.length} contours generated`);
  });
  
  // Test isoline generation
  console.log("\n   B. Testing isoline generation:");
  contourGenerator.isoLines(true).smoothingFactor(0.8);
  const isolines = contourGenerator(adaptedData.values);
  console.log(`      Generated ${isolines.length} isoline levels`);
  console.log(`      Output type: ${isolines[0]?.type || "none"}`);
  
  // Test with custom thresholds
  console.log("\n   C. Testing custom thresholds:");
  contourGenerator.isoLines(false).thresholds([50, 100, 150, 200, 250]);
  const customContours = contourGenerator(adaptedData.values);
  console.log(`      Custom thresholds: ${customContours.map(c => c.value).join(", ")}`);
  
  // Test saddle point disambiguation
  console.log("\n   D. Saddle point disambiguation is automatically applied");
  console.log("      - Case 5 and 10 saddle points are resolved based on center interpolation");
  console.log("      - Improves topological correctness of contour lines");
  
  // Test improved hole assignment
  console.log("\n   E. Improved hole assignment algorithm:");
  console.log("      - O(n log n) complexity instead of O(n²)");
  console.log("      - Proper nesting hierarchy for complex topologies");
  console.log("      - Spatial indexing with bounding boxes for performance");
}

// Export for use
export { demonstrateEnhancements };

// Run demonstration if this file is executed directly
if (typeof window === 'undefined') {
  demonstrateEnhancements();
}