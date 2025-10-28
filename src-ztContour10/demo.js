import contours from "./enhanced-contours.js";
import { contourData as testData1 } from "./testData1.js";
import { contourData as testData2 } from "./testData2.js";

/**
 * Enhanced D3-Contour Demo
 * 
 * This demo showcases the improved features of the enhanced d3-contour library:
 * 1. Saddle point disambiguation for better topology
 * 2. Optimized hole assignment algorithm (O(n) vs O(n²))
 * 3. Enhanced null value handling
 * 4. Adjustable smoothing parameters
 * 5. GeoJSON export functionality
 */
// Simple GeoJSON export function
function exportGeoJSON(contourLines) {
  const featureCollection = {
    type: "FeatureCollection",
    features: contourLines.map(contour => ({
      type: "Feature",
      properties: contour.properties || { value: contour.value },
      geometry: contour
    }))
  };
  return featureCollection;
}
function runDemo() {
  console.log("=== Enhanced D3-Contour Demo ===\n");

  // Demo 1: Process testData1 with null values
  console.log("Demo 1: Processing testData1 (contains null values)");
  demonstrateNullHandling();

  // Demo 2: Process testData2 with enhanced features
  console.log("\nDemo 2: Processing testData2 with enhanced features");
  demonstrateEnhancedFeatures();

  // Demo 3: Performance comparison
  console.log("\nDemo 3: Performance comparison");
  demonstratePerformance();

  console.log("\n=== Demo Complete ===");
}

function demonstrateNullHandling() {
  const { x, y, v } = testData1.data;
  const dx = x.length;
  const dy = y.length;
  
  // Flatten 2D array
  const values = [];
  for (let j = 0; j < dy; j++) {
    for (let i = 0; i < dx; i++) {
      values.push(v[j][i]);
    }
  }

  // Count null values
  const nullCount = values.filter(val => val == null).length;
  const totalCount = values.length;
  
  console.log(`  Dataset: ${dx} x ${dy} grid (${totalCount} cells)`);
  console.log(`  Null values: ${nullCount} (${(nullCount/totalCount*100).toFixed(1)}%)`);

  // Create enhanced contour generator
  const contourGen = contours()
    .size([dx, dy])
    .nullHandling(true)
    .smooth(true)
    .smoothFactor(0.7)
    .saddleDisambiguation(true)
    .geoJSON(true);

  // Generate contours with automatic thresholds
  const contourLines = contourGen(values);
  
  console.log(`  Generated: ${contourLines.length} contour levels`);
  
  // Show contour statistics
  contourLines.forEach((contour, i) => {
    const polygonCount = contour.coordinates.length;
    const totalRings = contour.coordinates.reduce((sum, poly) => sum + poly.length, 0);
    console.log(`    Level ${contour.value.toFixed(1)}: ${polygonCount} polygons, ${totalRings} rings`);
  });
}

function demonstrateEnhancedFeatures() {
  const { x, y, v } = testData2.data;
  const dx = x.length;
  const dy = y.length;
  
  // Flatten 2D array
  const values = [];
  for (let j = 0; j < dy; j++) {
    for (let i = 0; i < dx; i++) {
      values.push(v[j][i]);
    }
  }


       [
      [0, 'rgba(126,84,255,255)'],
      [10, 'rgba(99,66,255,255)'],
      [20, 'rgba(72,48,255,255)'],
      [30, 'rgba(45,30,255,255)'],
      [40, 'rgba(18,12,255,255)'],
      [50, 'rgba(0,15,240,255)'],
      [60, 'rgba(0,60,195,255)'],
      [70, 'rgba(0,105,150,255)'],
      [80, 'rgba(0,150,105,255)'],
      [90, 'rgba(0,195,60,255)'],
      [100, 'rgba(0,240,15,255)'],
      [110, 'rgba(30,255,0,255)'],
      [120, 'rgba(75,255,0,255)'],
      [130, 'rgba(120,255,0,255)'],
      [140, 'rgba(165,255,0,255)'],
      [150, 'rgba(210,255,0,255)'],
      [160, 'rgba(255,255,0,255)'],
      [170, 'rgba(255,228,0,255)'],
      [180, 'rgba(255,201,0,255)'],
      [190, 'rgba(255,174,0,255)'],
      [200, 'rgba(255,147,0,255)'],
      [210, 'rgba(255,120,0,255)'],
      [220, 'rgba(255,96,0,255)'],
      [230, 'rgba(255,78,0,255)'],
      [240, 'rgba(255,60,0,255)'],
      [250, 'rgba(255,42,0,255)'],
      [260, 'rgba(255,24,0,255)'],
      [270, 'rgba(255,6,0,255)'],
    ]

  console.log(`  Dataset: ${dx} x ${dy} grid (${values.length} cells)`);

  // Create contour generator with all enhanced features
  const contourGen = contours()
    .size([dx, dy])
    .smooth(true)
    .smoothFactor(1.8)
    .saddleDisambiguation(true)
    .geoJSON(true).thresholds([50, 100, 150, 200, 250, 300, 350]);

  const start = performance.now();
  const contourLines = contourGen(values);
  const duration = performance.now() - start;

  const jj = exportGeoJSON(contourLines);
  console.log(jj);
  console.log(`  Processing time: ${duration.toFixed(2)}ms`);
  console.log(`  Generated: ${contourLines.length} contour levels`);

  // Analyze contour quality
  let totalPolygons = 0;
  let totalHoles = 0;
  
  contourLines.forEach(contour => {
    totalPolygons += contour.coordinates.length;
    contour.coordinates.forEach(polygon => {
      totalHoles += Math.max(0, polygon.length - 1);
    });
  });

  console.log(`  Total polygons: ${totalPolygons}`);
  console.log(`  Total holes: ${totalHoles}`);
  console.log(`  Avg polygons per level: ${(totalPolygons / contourLines.length).toFixed(1)}`);

  // Demonstrate GeoJSON export
  const geoJsonData = exportGeoJSON(contourLines);
  
  console.log(`  GeoJSON export: ${(geoJsonData.length / 1024).toFixed(1)}KB`);
}

function demonstratePerformance() {
  // Create test datasets of different sizes
  const sizes = [20, 50, 100];
  
  sizes.forEach(size => {
    console.log(`\n  Testing ${size}x${size} grid:`);
    
    // Generate test data
    const values = new Array(size * size);
    for (let i = 0; i < values.length; i++) {
      const x = (i % size) / size * 10;
      const y = Math.floor(i / size) / size * 10;
      values[i] = Math.sin(x) * Math.cos(y) * 100 + Math.random() * 20;
    }

    // Test with basic settings
    const basicGen = contours()
      .size([size, size])
      .smooth(false)
      .saddleDisambiguation(false);

    const basicStart = performance.now();
    const basicResult = basicGen.thresholds(5)(values);
    const basicDuration = performance.now() - basicStart;

    // Test with enhanced settings
    const enhancedGen = contours()
      .size([size, size])
      .smooth(true)
      .smoothFactor(0.5)
      .saddleDisambiguation(true);

    const enhancedStart = performance.now();
    const enhancedResult = enhancedGen.thresholds(5)(values);
    const enhancedDuration = performance.now() - enhancedStart;

    console.log(`    Basic mode: ${basicDuration.toFixed(2)}ms, ${basicResult.length} levels`);
    console.log(`    Enhanced mode: ${enhancedDuration.toFixed(2)}ms, ${enhancedResult.length} levels`);
    console.log(`    Performance ratio: ${(enhancedDuration / basicDuration).toFixed(2)}x`);
  });
}

// Configuration examples
export const presets = {
  // Fast processing for real-time applications
  realtime: {
    smooth: false,
    smoothFactor: 0,
    saddleDisambiguation: false,
    geoJSON: false
  },
  
  // Balanced quality and performance
  balanced: {
    smooth: true,
    smoothFactor: 0.5,
    saddleDisambiguation: true,
    geoJSON: true
  },
  
  // High quality for final output
  highQuality: {
    smooth: true,
    smoothFactor: 0.8,
    saddleDisambiguation: true,
    geoJSON: true
  }
};

export function createContourGenerator(preset = 'balanced') {
  const config = presets[preset] || presets.balanced;
  
  return contours()
    .smooth(config.smooth)
    .smoothFactor(config.smoothFactor)
    .saddleDisambiguation(config.saddleDisambiguation)
    .geoJSON(config.geoJSON)
    .nullHandling(true);
}

// Run demo if executed directly
if (typeof window === 'undefined' && typeof global !== 'undefined') {
  runDemo();
}

export default runDemo;