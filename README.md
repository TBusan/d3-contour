# d3-contour (Enhanced)

<a href="https://d3js.org"><img src="https://github.com/d3/d3/raw/main/docs/public/logo.svg" width="256" height="256"></a>

This module computes contour polygons by applying [marching squares](https://en.wikipedia.org/wiki/Marching_squares) to a rectangular array of numeric values.

## 🆕 Enhanced Features

This enhanced version of d3-contour follows Plotly.js architecture to provide a complete contour visualization system with:

- **✅ Isoband Support** - Generate filled contour regions (not just lines)
- **✅ Geometric Consistency** - Isolines extracted from band boundaries (perfect alignment)
- **✅ Path Smoothing** - Chaikin and Catmull-Rom algorithms for smooth curves
- **✅ Label Engine** - Automatic label positioning, rotation, and collision detection
- **✅ Auto Levels** - Nice numbers algorithm for automatic threshold selection
- **✅ Dual Rendering** - Both SVG and Canvas renderers (with SSR support)
- **✅ Grid Padding** - Fixes boundary artifacts in contour generation

## Quick Start

```javascript
import { generateContoursAuto, renderContoursSVG } from 'd3-contour';

// Create grid data
const nx = 100, ny = 100;
const grid = new Float32Array(nx * ny);
for (let y = 0; y < ny; y++) {
  for (let x = 0; x < nx; x++) {
    grid[y * nx + x] = Math.sin(x / 10) * Math.cos(y / 10);
  }
}

// Generate contours with automatic levels
const { bands, lines } = generateContoursAuto(grid, {
  size: [nx, ny],
  count: 10
});

// Render to SVG
const svg = renderContoursSVG({ bands, lines }, {
  width: 800,
  height: 600
});
```

## Documentation

- **[API Documentation](./docs/API.md)** - Complete API reference
- **[Examples](./examples/)** - Usage examples
- **[Original d3-contour docs](https://d3js.org/d3-contour)** - Base API

## Architecture

```
Grid Data
  ↓
d3-contour          ← Computes isobands (contour regions)
  ↓
Extract Boundaries  ← Derives isolines from band edges
  ↓
Path Smoothing      ← Geometric smoothing (Chaikin, Catmull-Rom)
  ↓
Label Engine        ← Positioning, rotation, collision detection
  ↓
Renderer            ← SVG or Canvas output
```

### Key Principle

**Isolines are NOT computed independently.** They are extracted from isoband boundaries, ensuring perfect alignment between contour lines and filled regions - just like Plotly.js.

## Installation

```bash
npm install d3-contour
```

## Features Overview

### 1. Isoband Generation

Generate filled contour regions with proper topology:

```javascript
import { generateContours } from 'd3-contour';

const { bands, lines } = generateContours(grid, {
  size: [100, 100],
  thresholds: [0.2, 0.4, 0.6, 0.8]
});
```

### 2. Path Smoothing

Apply geometric smoothing for publication-quality curves:

```javascript
import { chaikinSmooth } from 'd3-contour';

const smoothed = lines.map(line => ({
  ...line,
  paths: line.paths.map(path => chaikinSmooth(path, 3))
}));
```

### 3. Contour Labels

Automatic label placement with rotation and collision detection:

```javascript
import {
  generateLabelsForContours,
  removeCollidingLabels
} from 'd3-contour';

let labels = generateLabelsForContours(lines);
labels = removeCollidingLabels(labels, {
  width: 60,
  height: 20,
  padding: 10
});
```

### 4. Auto Levels

Generate nice threshold values automatically:

```javascript
import { autoLevels } from 'd3-contour';

const levels = autoLevels(min, max, { count: 10 });
// Returns: [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
```

### 5. Dual Rendering

Render to SVG or Canvas with the same API:

```javascript
import { renderContoursSVG, renderContoursCanvas } from 'd3-contour';

// SVG
const svg = renderContoursSVG({ bands, lines }, options);

// Canvas
const canvas = renderToCanvas({ bands, lines }, 800, 600, options);
```

## Resources

- [Documentation](https://d3js.org/d3-contour)
- [Examples](https://observablehq.com/collection/@d3/d3-contour)
- [Releases](https://github.com/d3/d3-contour/releases)
- [Getting help](https://d3js.org/community)
