/**
 * Basic usage example for enhanced d3-contour
 *
 * This example demonstrates the core workflow for generating
 * contours with proper padding, smoothing, and rendering.
 */

import {
  // Data generation
  generateContours,
  generateContoursAuto,
  padGrid,

  // Smoothing
  chaikinSmooth,

  // Level generation
  autoLevels,

  // Labels
  generateLabelsForContours,
  removeCollidingLabels,

  // Rendering
  renderContoursSVG,
  renderContoursCanvas,
  createSVGColorScale
} from '../src/index.js';

// Example 1: Basic contour generation with manual levels
// =======================================================

function example1_manualLevels() {
  // Create sample grid data (100x100)
  const nx = 100;
  const ny = 100;
  const grid = new Float32Array(nx * ny);

  // Generate some interesting data (e.g., a mixture of Gaussians)
  for (let y = 0; y < ny; y++) {
    for (let x = 0; x < nx; x++) {
      const value =
        Math.exp(-((x - 30) ** 2 + (y - 30) ** 2) / 400) +
        Math.exp(-((x - 70) ** 2 + (y - 70) ** 2) / 400) +
        Math.sin(x / 10) * Math.cos(y / 10);
      grid[y * nx + x] = value;
    }
  }

  // Option 1: Without padding (will have boundary artifacts)
  const { bands, lines } = generateContours(grid, {
    size: [nx, ny],
    thresholds: [0.2, 0.4, 0.6, 0.8, 1.0]
  });

  // Option 2: With padding (recommended)
  const paddedGrid = padGrid(grid, nx, ny);
  const { bands: paddedBands, lines: paddedLines } = generateContours(paddedGrid, {
    size: [nx + 2, ny + 2],
    thresholds: [0.2, 0.4, 0.6, 0.8, 1.0]
  });

  console.log(`Generated ${bands.length} bands and ${lines.length} lines`);
  return { bands, lines };
}

// Example 2: Automatic level selection
// ======================================

function example2_autoLevels() {
  const nx = 100;
  const ny = 100;
  const grid = new Float32Array(nx * ny);

  // Generate data
  for (let y = 0; y < ny; y++) {
    for (let x = 0; x < nx; x++) {
      grid[y * nx + x] = Math.sin(x / 10) * Math.cos(y / 10);
    }
  }

  // Automatic level generation
  const { bands, lines, thresholds } = generateContoursAuto(grid, {
    size: [nx, ny],
    count: 8
  });

  console.log('Auto-generated thresholds:', thresholds);
  return { bands, lines, thresholds };
}

// Example 3: Path smoothing
// ==========================

function example3_smoothing() {
  const nx = 100;
  const ny = 100;
  const grid = new Float32Array(nx * ny);

  for (let y = 0; y < ny; y++) {
    for (let x = 0; x < nx; x++) {
      grid[y * nx + x] = Math.sin(x / 5) * Math.cos(y / 5);
    }
  }

  const { bands, lines } = generateContours(grid, {
    size: [nx, ny],
    thresholds: [0.2, 0.4, 0.6, 0.8]
  });

  // Apply smoothing to all line paths
  const smoothedLines = lines.map(line => ({
    ...line,
    paths: line.paths.map(path => chaikinSmooth(path, 2))
  }));

  return { bands, lines: smoothedLines };
}

// Example 4: SVG rendering with labels
// =====================================

function example4_svgRendering() {
  const nx = 100;
  const ny = 100;
  const grid = new Float32Array(nx * ny);

  for (let y = 0; y < ny; y++) {
    for (let x = 0; x < nx; x++) {
      grid[y * nx + x] = Math.exp(-((x - 50) ** 2 + (y - 50) ** 2) / 500);
    }
  }

  const { bands, lines } = generateContours(grid, {
    size: [nx, ny],
    thresholds: [0.1, 0.3, 0.5, 0.7, 0.9]
  });

  // Generate labels
  let labels = generateLabelsForContours(lines);
  labels = removeCollidingLabels(labels, { width: 50, height: 15 });

  // Create color scale
  const colorScale = createSVGColorScale('viridis', [0, 1]);

  // Render to SVG
  const svg = renderContoursSVG(
    { bands, lines },
    {
      width: 800,
      height: 600,
      bandOptions: {
        fill: (value) => colorScale((value[0] + value[1]) / 2),
        opacity: 0.7
      },
      lineOptions: {
        stroke: (level) => level > 0.5 ? 'red' : 'black',
        strokeWidth: 1.5
      },
      labelOptions: {
        fontSize: 14,
        halo: 'white',
        haloWidth: 4
      },
      labels
    }
  );

  return svg;
}

// Example 5: Canvas rendering (browser)
// ======================================

function example5_canvasRendering() {
  const nx = 100;
  const ny = 100;
  const grid = new Float32Array(nx * ny);

  for (let y = 0; y < ny; y++) {
    for (let x = 0; x < nx; x++) {
      grid[y * nx + x] = Math.sin(x / 8) * Math.cos(y / 8);
    }
  }

  const { bands, lines } = generateContoursAuto(grid, {
    size: [nx, ny],
    count: 6
  });

  // Create canvas (browser)
  const canvas = document.createElement('canvas');
  canvas.width = 800;
  canvas.height = 600;
  document.body.appendChild(canvas);

  const ctx = canvas.getContext('2d');

  // Render
  const colorScale = createSVGColorScale('plasma', [-1, 1]);

  renderContoursCanvas(
    ctx,
    { bands, lines },
    {
      bandOptions: {
        fill: (value) => colorScale((value[0] + value[1]) / 2),
        opacity: 0.6
      },
      lineOptions: {
        stroke: () => 'rgba(0,0,0,0.3)',
        strokeWidth: 1
      }
    }
  );

  return canvas;
}

// Example 6: Complete workflow with all features
// ===============================================

function example6_completeWorkflow() {
  const nx = 100;
  const ny = 100;
  const grid = new Float32Array(nx * ny);

  // Generate test data
  for (let y = 0; y < ny; y++) {
    for (let x = 0; x < nx; x++) {
      const dx = x - 50;
      const dy = y - 50;
      grid[y * nx + x] = Math.exp(-(dx * dx + dy * dy) / 400);
    }
  }

  // Step 1: Pad grid to fix boundary artifacts
  const paddedGrid = padGrid(grid, nx, ny);

  // Step 2: Generate contours with automatic levels
  const { bands, lines, thresholds } = generateContoursAuto(paddedGrid, {
    size: [nx + 2, ny + 2],
    count: 10
  });

  // Step 3: Apply smoothing to lines
  const smoothedLines = lines.map(line => ({
    ...line,
    paths: line.paths.map(path => chaikinSmooth(path, 3))
  }));

  // Step 4: Generate labels and remove collisions
  const labels = removeCollidingLabels(
    generateLabelsForContours(smoothedLines),
    { width: 60, height: 20, padding: 10 }
  );

  // Step 5: Render to SVG
  const colorScale = createSVGColorScale('viridis', [0, 1]);

  const svg = renderContoursSVG(
    { bands, lines: smoothedLines },
    {
      width: 800,
      height: 600,
      bandOptions: {
        fill: (value) => colorScale((value[0] + value[1]) / 2),
        opacity: 0.7
      },
      lineOptions: {
        stroke: () => 'rgba(255,255,255,0.8)',
        strokeWidth: 2
      },
      labelOptions: {
        fontSize: 12,
        fill: 'black',
        halo: 'white',
        haloWidth: 3
      },
      labels
    }
  );

  return { svg, thresholds };
}

// Example 7: Server-side rendering (Node.js)
// ===========================================

async function example7_ssr() {
  const nx = 100;
  const ny = 100;
  const grid = new Float32Array(nx * ny);

  for (let y = 0; y < ny; y++) {
    for (let x = 0; x < nx; x++) {
      grid[y * nx + x] = Math.sin(x / 10) * Math.cos(y / 10);
    }
  }

  const { bands, lines } = generateContoursAuto(grid, {
    size: [nx, ny],
    count: 8
  });

  const colorScale = createSVGColorScale('blues', [-1, 1]);

  const svg = renderContoursSVG(
    { bands, lines },
    {
      width: 1200,
      height: 900,
      bandOptions: {
        fill: (value) => colorScale((value[0] + value[1]) / 2),
        opacity: 0.8
      },
      lineOptions: {
        stroke: () => 'white',
        strokeWidth: 1.5
      }
    }
  );

  // Save SVG to file
  const fs = await import('fs/promises');
  await fs.writeFile('contour-output.svg', svg);

  console.log('SVG saved to contour-output.svg');

  // For PNG output (requires @napi-rs/canvas):
  // import { createCanvas } from '@napi-rs/canvas';
  // import { renderContoursToPNG } from '../src/index.js';
  // await renderContoursToPNG(
  //   { bands, lines },
  //   'contour-output.png',
  //   1200,
  //   900,
  //   { bandOptions: { fill: colorScale, opacity: 0.8 } }
  // );

  return svg;
}

// Example 8: Using nice levels for better visualization
// =====================================================

function example8_niceLevels() {
  const nx = 100;
  const ny = 100;
  const grid = new Float32Array(nx * ny);

  for (let y = 0; y < ny; y++) {
    for (let x = 0; x < nx; x++) {
      grid[y * nx + x] = (x + y) / 200 + Math.random() * 0.1;
    }
  }

  // Calculate data range
  const min = Math.min(...grid);
  const max = Math.max(...grid);

  // Generate nice levels
  const levels = autoLevels(min, max, { count: 10 });

  console.log('Nice levels:', levels);

  const { bands, lines } = generateContours(grid, {
    size: [nx, ny],
    thresholds: levels
  });

  return { bands, lines, levels };
}

// Export examples
export {
  example1_manualLevels,
  example2_autoLevels,
  example3_smoothing,
  example4_svgRendering,
  example5_canvasRendering,
  example6_completeWorkflow,
  example7_ssr,
  example8_niceLevels
};

// Run examples if this file is executed directly
if (import.meta.url === `file://${process.argv[1]}`) {
  console.log('Running examples...');

  console.log('\n=== Example 1: Manual Levels ===');
  example1_manualLevels();

  console.log('\n=== Example 2: Auto Levels ===');
  example2_autoLevels();

  console.log('\n=== Example 8: Nice Levels ===');
  example8_niceLevels();

  console.log('\nExamples complete!');
}
