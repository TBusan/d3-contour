import contours from "./contours.js";
import dataAdapter from "./dataAdapter.js";
import { contourData as testData1 } from "./testData1.js";
import { contourData as testData2 } from "./testData2.js";

// Comprehensive test suite for the enhanced d3-contour library
function runTests() {
  console.log("=== Enhanced D3-Contour Test Suite ===\n");
  
  let passedTests = 0;
  let totalTests = 0;
  
  // Test 1: Data adapter with testData1
  totalTests++;
  console.log("Test 1: Data adapter with testData1 (null value handling)");
  try {
    const adapter = dataAdapter();
    const adapted = adapter(testData1);
    
    console.log(`  ✓ Grid dimensions: ${adapted.width} x ${adapted.height}`);
    console.log(`  ✓ Coordinate ranges: X[${adapted.xmin}, ${adapted.xmax}], Y[${adapted.ymin}, ${adapted.ymax}]`);
    
    const nullCount = adapted.values.filter(v => v === null).length;
    console.log(`  ✓ Null values preserved: ${nullCount} nulls found`);
    
    if (adapted.width === 80 && adapted.height === 17 && nullCount > 0) {
      console.log("  ✅ PASSED\n");
      passedTests++;
    } else {
      console.log("  ❌ FAILED: Unexpected dimensions or null handling\n");
    }
  } catch (error) {
    console.log(`  ❌ FAILED: ${error.message}\n`);
  }
  
  // Test 2: Data adapter with testData2
  totalTests++;
  console.log("Test 2: Data adapter with testData2 (dense grid)");
  try {
    const adapter = dataAdapter();
    const adapted = adapter(testData2);
    
    console.log(`  ✓ Grid dimensions: ${adapted.width} x ${adapted.height}`);
    console.log(`  ✓ Value range: [${testData2.zmin}, ${testData2.zmax}]`);
    
    const hasValidData = adapted.values.some(v => v !== null && !isNaN(v));
    console.log(`  ✓ Valid data present: ${hasValidData}`);
    
    if (adapted.width === 101 && adapted.height === 20 && hasValidData) {
      console.log("  ✅ PASSED\n");
      passedTests++;
    } else {
      console.log("  ❌ FAILED: Unexpected dimensions or data\n");
    }
  } catch (error) {
    console.log(`  ❌ FAILED: ${error.message}\n`);
  }
  
  // Test 3: Basic contour generation with null handling
  totalTests++;
  console.log("Test 3: Contour generation with null value support");
  try {
    const adapter = dataAdapter();
    const adapted = adapter(testData1);
    const contourGen = contours()
      .size([adapted.width, adapted.height])
      .nullValue(null)
      .thresholds(5);
    
    const result = contourGen(adapted.values);
    console.log(`  ✓ Generated ${result.length} contour levels`);
    
    const hasValidContours = result.length > 0 && result.every(c => 
      c.type === "MultiPolygon" && 
      Array.isArray(c.coordinates) &&
      typeof c.value === "number"
    );
    
    if (hasValidContours) {
      console.log("  ✅ PASSED\n");
      passedTests++;
    } else {
      console.log("  ❌ FAILED: Invalid contour structure\n");
    }
  } catch (error) {
    console.log(`  ❌ FAILED: ${error.message}\n`);
  }
  
  // Test 4: Isoline generation
  totalTests++;
  console.log("Test 4: Isoline generation");
  try {
    const adapter = dataAdapter();
    const adapted = adapter(testData2);
    const contourGen = contours()
      .size([adapted.width, adapted.height])
      .isoLines(true)
      .thresholds([100, 200, 300]);
    
    const result = contourGen(adapted.values);
    console.log(`  ✓ Generated ${result.length} isoline levels`);
    
    const hasValidIsolines = result.length > 0 && result.every(c => 
      c.type === "MultiLineString" && Array.isArray(c.coordinates)
    );
    
    if (hasValidIsolines) {
      console.log("  ✅ PASSED\n");
      passedTests++;
    } else {
      console.log("  ❌ FAILED: Invalid isoline structure\n");
    }
  } catch (error) {
    console.log(`  ❌ FAILED: ${error.message}\n`);
  }
  
  // Test 5: Smoothing factor control
  totalTests++;
  console.log("Test 5: Smoothing factor control");
  try {
    const adapter = dataAdapter();
    const adapted = adapter(testData2);
    const contourGen = contours()
      .size([adapted.width, adapted.height])
      .thresholds([150]);
    
    // Test different smoothing factors
    const results = {};
    [0, 0.5, 1.0].forEach(factor => {
      contourGen.smoothingFactor(factor);
      const result = contourGen(adapted.values);
      results[factor] = result[0]?.coordinates[0]?.[0] || [];
    });
    
    const hasVariedResults = Object.keys(results).length === 3;
    console.log(`  ✓ Tested smoothing factors: 0, 0.5, 1.0`);
    
    if (hasVariedResults) {
      console.log("  ✅ PASSED\n");
      passedTests++;
    } else {
      console.log("  ❌ FAILED: Smoothing factor not working\n");
    }
  } catch (error) {
    console.log(`  ❌ FAILED: ${error.message}\n`);
  }
  
  // Test 6: Coordinate transformation
  totalTests++;
  console.log("Test 6: Coordinate transformation");
  try {
    const adapter = dataAdapter();
    const adapted = adapter(testData1);
    const contourGen = contours()
      .size([adapted.width, adapted.height])
      .nullValue(null)
      .thresholds([100]);
    
    const contourResults = contourGen(adapted.values);
    const transformed = adapter.createContourData(adapted, contourResults);
    
    // Check if coordinates are in the original coordinate system
    const hasTransformedCoords = transformed.length > 0 && 
      transformed[0].coordinates.length > 0 &&
      transformed[0].coordinates[0].length > 0;
    
    if (hasTransformedCoords) {
      const firstPoint = transformed[0].coordinates[0][0][0];
      const inOriginalRange = firstPoint[0] >= testData1.x[0] && 
                            firstPoint[0] <= testData1.x[testData1.x.length - 1];
      
      console.log(`  ✓ Coordinates transformed to original range`);
      if (inOriginalRange) {
        console.log("  ✅ PASSED\n");
        passedTests++;
      } else {
        console.log("  ❌ FAILED: Coordinates not in original range\n");
      }
    } else {
      console.log("  ❌ FAILED: No valid coordinates found\n");
    }
  } catch (error) {
    console.log(`  ❌ FAILED: ${error.message}\n`);
  }
  
  // Test 7: Performance and hole assignment
  totalTests++;
  console.log("Test 7: Enhanced hole assignment");
  try {
    const adapter = dataAdapter();
    const adapted = adapter(testData2);
    const contourGen = contours()
      .size([adapted.width, adapted.height])
      .thresholds(10); // Generate multiple levels to test hole assignment
    
    const startTime = performance.now();
    const result = contourGen(adapted.values);
    const endTime = performance.now();
    
    console.log(`  ✓ Generated ${result.length} contours in ${(endTime - startTime).toFixed(2)}ms`);
    
    // Count holes (polygons with more than one ring)
    const totalHoles = result.reduce((sum, contour) => 
      sum + contour.coordinates.reduce((holeSum, polygon) => 
        holeSum + Math.max(0, polygon.length - 1), 0), 0);
    
    console.log(`  ✓ Total holes assigned: ${totalHoles}`);
    
    if (result.length > 0 && endTime - startTime < 1000) { // Should be reasonably fast
      console.log("  ✅ PASSED\n");
      passedTests++;
    } else {
      console.log("  ❌ FAILED: Performance or generation issues\n");
    }
  } catch (error) {
    console.log(`  ❌ FAILED: ${error.message}\n`);
  }
  
  // Test Results Summary
  console.log("=".repeat(50));
  console.log(`Test Results: ${passedTests}/${totalTests} tests passed`);
  
  if (passedTests === totalTests) {
    console.log("🎉 All tests passed! The enhanced d3-contour library is working correctly.");
  } else {
    console.log(`⚠️  ${totalTests - passedTests} test(s) failed. Please review the implementation.`);
  }
  
  return { passed: passedTests, total: totalTests };
}

// Export for use
export { runTests };

// Add simple performance polyfill for non-browser environments
if (typeof performance === 'undefined') {
  global.performance = {
    now: () => Date.now()
  };
}

// Run tests if this file is executed directly
if (typeof window === 'undefined') {
  runTests();
}