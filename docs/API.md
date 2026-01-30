# Enhanced D3-Contour API Documentation

Complete API reference for the enhanced d3-contour library, following Plotly.js architecture principles.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Core Concepts](#core-concepts)
- [API Reference](#api-reference)
  - [Contour Generation](#contour-generation)
  - [Grid Padding](#grid-padding)
  - [Path Smoothing](#path-smoothing)
  - [Label Generation](#label-generation)
  - [Level Generation](#level-generation)
  - [SVG Rendering](#svg-rendering)
  - [Canvas Rendering](#canvas-rendering)
- [Examples](#examples)
- [Architecture](#architecture)

---

## Installation

```bash
npm install d3-contour
```

```javascript
// ES modules
import {
  generateContours,
  chaikinSmooth,
  renderContoursSVG
} from 'd3-contour';

// CommonJS
const {
  generateContours,
  chaikinSmooth,
  renderContoursSVG
} = require('d3-contour');
```

---

## Quick Start

```javascript
import { generateContoursAuto, renderContoursSVG } from 'd3-contour';

// 1. Create grid data
const nx = 100, ny = 100;
const grid = new Float32Array(nx * ny);
for (let y = 0; y < ny; y++) {
  for (let x = 0; x < nx; x++) {
    grid[y * nx + x] = Math.sin(x / 10) * Math.cos(y / 10);
  }
}

// 2. Generate contours with automatic levels
const { bands, lines } = generateContoursAuto(grid, {
  size: [nx, ny],
  count: 10
});

// 3. Render to SVG
const svg = renderContoursSVG({ bands, lines }, {
  width: 800,
  height: 600
});

document.body.innerHTML = svg;
```

---

## Core Concepts

### The Plotly Architecture

This library follows Plotly's approach to contour visualization:

1. **Grid Data** → Your input 2D data
2. **d3-contour** → Computes isobands (contour regions)
3. **Extract Boundaries** → Derives isolines from band edges
4. **Path Smoothing** → Applies geometric smoothing (not grid smoothing)
5. **Label Engine** → Positions and rotates labels
6. **Renderer** → Outputs SVG or Canvas

### Key Principle

**Isolines are NOT computed independently.** They are extracted from isoband boundaries, ensuring perfect alignment between lines and filled regions.

---

## API Reference

### Contour Generation

#### `generateContours(grid, options)`

Generate both isobands and isolines from grid data, ensuring perfect alignment.

**Parameters:**
- `grid` (Float32Array|Array) - Flat grid data [y * nx + x]
- `options` (Object)
  - `size` ([number, number]) - Grid dimensions [width, height]
  - `thresholds` (Array<number>) - Contour level values
  - `smooth` (boolean) - Enable built-in smoothing (default: false, use post-processing)

**Returns:**
```javascript
{
  bands: Array<{
    value: [number, number],  // [min, max] range
    polygons: Array,           // GeoJSON MultiPolygon format
    level: number              // Upper boundary level
  }>,
  lines: Array<{
    paths: Array<Array<[number, number]>>, // Line paths
    level: number                          // Level value
  }>
}
```

**Example:**
```javascript
const { bands, lines } = generateContours(grid, {
  size: [100, 100],
  thresholds: [0.2, 0.4, 0.6, 0.8]
});
```

---

#### `generateContoursAuto(grid, options)`

Generate contours with automatic level selection.

**Parameters:**
- `grid` (Float32Array|Array) - Flat grid data
- `options` (Object)
  - `size` ([number, number]) - Grid dimensions
  - `count` (number) - Target number of levels (default: 10)
  - `smooth` (boolean) - Enable smoothing (default: false)

**Returns:**
Same as `generateContours`, plus:
```javascript
{
  ...,
  thresholds: Array<number>  // Computed threshold values
}
```

**Example:**
```javascript
const { bands, lines, thresholds } = generateContoursAuto(grid, {
  size: [100, 100],
  count: 8
});
console.log('Auto levels:', thresholds);
```

---

#### `generateIsobands(grid, options)`

Generate only isobands (filled regions).

**Parameters:**
- `grid` (Float32Array|Array) - Flat grid data
- `options` (Object)
  - `size` ([number, number]) - Grid dimensions
  - `thresholds` (Array<[number, number]>) - Array of [min, max] ranges
  - `smooth` (boolean) - Enable smoothing

**Returns:**
```javascript
Array<{
  value: [number, number],
  polygons: Array,
  level: number
}>
```

**Example:**
```javascript
const bands = generateIsobands(grid, {
  size: [100, 100],
  thresholds: [[0, 0.5], [0.5, 1.0], [1.0, 1.5]]
});
```

---

#### `extractIsolinesFromBands(bands)`

Extract isolines from isoband boundaries.

**Parameters:**
- `bands` (Array) - Array of isoband objects

**Returns:**
```javascript
Array<{
  paths: Array<Array<[number, number]>>,
  level: number
}>
```

**Example:**
```javascript
const lines = extractIsolinesFromBands(bands);
```

---

### Grid Padding

#### `padGrid(grid, nx, ny)`

Pad grid to fix boundary artifacts in contour generation.

**Parameters:**
- `grid` (Float32Array|Array) - Original grid
- `nx` (number) - Grid width
- `ny` (number) - Grid height

**Returns:**
- `Float32Array` - Padded grid with size (nx + 2) * (ny + 2)

**Example:**
```javascript
const paddedGrid = padGrid(grid, 100, 100);
const contours = generateContours(paddedGrid, {
  size: [102, 102],
  thresholds: [0.5, 1.0]
});
```

---

#### `unpadGrid(padded, originalNx, originalNy)`

Remove padding from a padded grid.

**Parameters:**
- `padded` (Float32Array|Array) - Padded grid
- `originalNx` (number) - Original width
- `originalNy` (number) - Original height

**Returns:**
- `Float32Array` - Unpadded grid

---

#### `transformCoords(x, y, originalNx, originalNy)`

Transform coordinates from padded to original space.

**Parameters:**
- `x`, `y` (number) - Coordinates in padded space
- `originalNx`, `originalNy` (number) - Original dimensions

**Returns:**
- `[number, number]` - Coordinates in original space

---

### Path Smoothing

#### `chaikinSmooth(points, iterations)`

Apply Chaikin's corner cutting algorithm for path smoothing.

**Parameters:**
- `points` (Array<[number, number]>) - Path points
- `iterations` (number) - Smoothing iterations (default: 2)

**Returns:**
- `Array<[number, number]>` - Smoothed points

**Example:**
```javascript
const smoothedPath = chaikinSmooth(path, 3);
```

---

#### `catmullRomSmooth(points, tension, segments)`

Apply Catmull-Rom spline interpolation.

**Parameters:**
- `points` (Array<[number, number]>) - Path points
- `tension` (number) - Tension parameter (default: 0.5)
- `segments` (number) - Segments per curve (default: 10)

**Returns:**
- `Array<[number, number]>` - Smoothed points

**Example:**
```javascript
const smoothed = catmullRomSmooth(path, 0.5, 15);
```

---

#### `catmullRomToBezier(points, tension)`

Convert Catmull-Rom points to cubic Bezier control points for SVG.

**Parameters:**
- `points` (Array<[number, number]>) - Path points
- `tension` (number) - Tension parameter

**Returns:**
```javascript
Array<{
  p0: [number, number],
  cp1: [number, number],
  cp2: [number, number],
  p1: [number, number]
}>
```

---

#### `isPathClosed(points, tolerance)`

Check if a path is closed (polygon) or open (polyline).

**Parameters:**
- `points` (Array<[number, number]>) - Path points
- `tolerance` (number) - Distance tolerance (default: 1e-6)

**Returns:**
- `boolean`

---

#### `pathLength(points)`

Calculate the total length of a path.

**Parameters:**
- `points` (Array<[number, number]>) - Path points

**Returns:**
- `number` - Path length

---

#### `samplePathAt(points, distance)`

Sample a point at a specific distance along a path.

**Parameters:**
- `points` (Array<[number, number]>) - Path points
- `distance` (number) - Distance along path

**Returns:**
- `[number, number]` - Interpolated point

---

### Label Generation

#### `generateLabel(path, level, options)`

Generate a single label for a contour path.

**Parameters:**
- `path` (Array<[number, number]>) - Contour path
- `level` (number) - Contour level
- `options` (Object)
  - `format` (string) - Number format (default: '.2f')
  - `position` (number) - Position along path 0-1 (default: 0.5)

**Returns:**
```javascript
{
  x: number,
  y: number,
  angle: number,      // Rotation in radians
  text: string,
  level: number,
  pathLength: number
} | null
```

**Example:**
```javascript
const label = generateLabel(path, 1.5, {
  format: '.1f',
  position: 0.5  // Center of path
});
```

---

#### `generateLabels(path, level, options)`

Generate multiple labels for a single contour path.

**Parameters:**
- `path` (Array<[number, number]>) - Contour path
- `level` (number) - Contour level
- `options` (Object)
  - `maxLabels` (number) - Maximum labels (default: 3)
  - `minSpacing` (number) - Minimum spacing (default: 50)
  - `format` (string) - Number format

**Returns:**
- `Array<Object>` - Array of label objects

---

#### `generateLabelsForContours(contours, options)`

Generate labels for multiple contours.

**Parameters:**
- `contours` (Array) - Array of `{paths, level}` objects
- `options` (Object) - Label options

**Returns:**
- `Array<Object>` - All label objects

**Example:**
```javascript
const labels = generateLabelsForContours(lines);
```

---

#### `removeCollidingLabels(labels, options)`

Remove colliding labels using greedy algorithm.

**Parameters:**
- `labels` (Array) - Array of label objects
- `options` (Object)
  - `width` (number) - Label width in pixels (default: 60)
  - `height` (number) - Label height (default: 20)
  - `padding` (number) - Padding around labels (default: 5)

**Returns:**
- `Array<Object>` - Filtered labels without collisions

**Example:**
```javascript
const cleanLabels = removeCollidingLabels(labels, {
  width: 50,
  height: 15,
  padding: 10
});
```

---

#### `removeCollidingLabelsSpatial(labels, options)`

Remove colliding labels using spatial indexing (faster for many labels).

**Parameters:**
- Same as `removeCollidingLabels`, plus:
  - `cellSize` (number) - Spatial grid cell size (default: 70)

**Returns:**
- `Array<Object>` - Filtered labels

---

#### `formatLevel(level, format)`

Format a level value as a string.

**Parameters:**
- `level` (number) - Level value
- `format` (string) - D3 format specifier (e.g., '.2f', '.1s')

**Returns:**
- `string` - Formatted string

---

### Level Generation

#### `autoLevels(min, max, options)`

Generate nice levels for a data range.

**Parameters:**
- `min`, `max` (number) - Data range
- `options` (Object)
  - `count` (number) - Target number of levels (default: 10)
  - `includeEndpoints` (boolean) - Include min/max (default: false)
  - `minLevels` (number) - Minimum number of levels (default: 2)

**Returns:**
- `Array<number>` - Nice level values

**Example:**
```javascript
const levels = autoLevels(3.7, 97.3, { count: 8 });
// Returns: [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
```

---

#### `niceStep(range, targetSteps)`

Find a nice step size for a range.

**Parameters:**
- `range` (number) - Data range (max - min)
- `targetSteps` (number) - Target number of steps

**Returns:**
- `number` - Nice step size

---

#### `quantileLevels(data, options)`

Generate levels using quantiles (for uneven data distribution).

**Parameters:**
- `data` (Array|Float32Array) - Data array
- `options` (Object)
  - `count` (number) - Target number of levels

**Returns:**
- `Array<number>` - Level values at quantiles

---

#### `logLevels(min, max, options)`

Generate logarithmic levels for exponential data.

**Parameters:**
- `min`, `max` (number) - Range (must be positive)
- `options` (Object)
  - `count` (number) - Number of levels

**Returns:**
- `Array<number>` - Logarithmically spaced levels

---

#### `levelsToBands(levels, options)`

Convert level values to band ranges.

**Parameters:**
- `levels` (Array<number>) - Level values
- `options` (Object)
  - `min` (number) - Minimum data value
  - `max` (number) - Maximum data value
  - `extend` (boolean) - Extend beyond min/max (default: true)

**Returns:**
- `Array<[number, number]>` - Band ranges

---

### SVG Rendering

#### `renderContoursSVG(contours, options)`

Render complete contour visualization as SVG.

**Parameters:**
- `contours` (Object) - `{bands, lines}`
- `options` (Object)
  - `width`, `height` (number) - SVG dimensions
  - `viewBox` (string) - SVG viewBox
  - `labels` (Array) - Optional labels
  - `bandOptions` (Object) - Band rendering options
    - `fill` (Function) - Color scale function
    - `opacity` (number) - Fill opacity
  - `lineOptions` (Object) - Line rendering options
    - `stroke` (Function) - Stroke color function
    - `strokeWidth` (number) - Line width
  - `labelOptions` (Object) - Label rendering options
    - `fontSize` (number) - Font size
    - `fill` (string) - Text color
    - `halo` (string) - Halo color

**Returns:**
- `string` - Complete SVG document

**Example:**
```javascript
const svg = renderContoursSVG({ bands, lines }, {
  width: 800,
  height: 600,
  bandOptions: {
    fill: (v) => `rgba(0,0,255,${v[0]/10})`,
    opacity: 0.7
  },
  lineOptions: {
    stroke: (level) => level > 0.5 ? 'red' : 'black',
    strokeWidth: 1.5
  }
});
```

---

#### `renderBandsSVG(bands, options)`

Render only isobands as SVG paths.

**Parameters:**
- `bands` (Array) - Isoband objects
- `options` (Object)
  - `fill` (Function) - Color function
  - `opacity` (number) - Fill opacity
  - `class` (string) - CSS class

**Returns:**
- `string` - SVG path elements

---

#### `renderLinesSVG(lines, options)`

Render only isolines as SVG paths.

**Parameters:**
- `lines` (Array) - Isoline objects
- `options` (Object)
  - `stroke` (Function) - Color function
  - `strokeWidth` (number) - Line width
  - `strokeLinecap` (string) - Line cap style
  - `strokeLinejoin` (string) - Line join style

**Returns:**
- `string` - SVG path elements

---

#### `renderLabelsSVG(labels, options)`

Render labels as SVG text elements.

**Parameters:**
- `labels` (Array) - Label objects
- `options` (Object)
  - `fontSize` (number)
  - `fontFamily` (string)
  - `fill` (string) - Text color
  - `halo` (string) - Halo color
  - `haloWidth` (number) - Halo width

**Returns:**
- `string` - SVG text elements

---

#### `createSVGColorScale(scheme, domain)`

Create a color scale function.

**Parameters:**
- `scheme` (string) - 'viridis', 'plasma', 'blues', 'heatmap'
- `domain` ([number, number]) - Value range

**Returns:**
- `Function` - Color scale function (value) => color string

**Example:**
```javascript
const colorScale = createSVGColorScale('viridis', [0, 1]);
const color = colorScale(0.5); // Returns 'rgb(...)'
```

---

### Canvas Rendering

#### `renderContoursCanvas(ctx, contours, options)`

Render contours on Canvas context.

**Parameters:**
- `ctx` (CanvasRenderingContext2D) - Canvas 2D context
- `contours` (Object) - `{bands, lines}`
- `options` (Object) - Similar to SVG options

**Example:**
```javascript
const canvas = document.createElement('canvas');
const ctx = canvas.getContext('2d');

renderContoursCanvas(ctx, { bands, lines }, {
  bandOptions: { fill: colorScale, opacity: 0.7 },
  lineOptions: { stroke: 'black', strokeWidth: 1 }
});
```

---

#### `renderBandsCanvas(ctx, bands, options)`

Render only isobands on Canvas.

---

#### `renderLinesCanvas(ctx, lines, options)`

Render only isolines on Canvas.

---

#### `renderLabelsCanvas(ctx, labels, options)`

Render labels on Canvas.

---

#### `renderToCanvas(contours, width, height, options)`

Create a new canvas and render contours on it.

**Parameters:**
- `contours` (Object) - `{bands, lines}`
- `width`, `height` (number) - Canvas dimensions
- `options` (Object) - Rendering options

**Returns:**
- `HTMLCanvasElement` - Canvas with rendered contours

---

#### `renderToDataURL(contours, width, height, options)`

Render contours and export as PNG data URL (browser only).

**Returns:**
- `string` - PNG data URL

---

#### `renderContoursToPNG(contours, filepath, width, height, options)`

Render contours and save as PNG (Node.js with @napi-rs/canvas).

**Example:**
```javascript
import { createCanvas } from '@napi-rs/canvas';

await renderContoursToPNG(
  { bands, lines },
  'output.png',
  800,
  600,
  { bandOptions: { opacity: 0.7 } }
);
```

---

## Examples

### Basic Usage

```javascript
import {
  generateContoursAuto,
  renderContoursSVG,
  createSVGColorScale
} from 'd3-contour';

// Create grid
const nx = 100, ny = 100;
const grid = new Float32Array(nx * ny);
for (let y = 0; y < ny; y++) {
  for (let x = 0; x < nx; x++) {
    grid[y * nx + x] = Math.sin(x / 10) * Math.cos(y / 10);
  }
}

// Generate contours
const { bands, lines } = generateContoursAuto(grid, {
  size: [nx, ny],
  count: 10
});

// Create color scale
const colorScale = createSVGColorScale('viridis', [-1, 1]);

// Render
const svg = renderContoursSVG({ bands, lines }, {
  width: 800,
  height: 600,
  bandOptions: {
    fill: (value) => colorScale((value[0] + value[1]) / 2),
    opacity: 0.7
  }
});

document.body.innerHTML = svg;
```

### With Smoothing and Labels

```javascript
import {
  generateContours,
  chaikinSmooth,
  generateLabelsForContours,
  removeCollidingLabels,
  renderContoursSVG
} from 'd3-contour';

// Generate contours
const { bands, lines } = generateContours(grid, {
  size: [100, 100],
  thresholds: [0.2, 0.4, 0.6, 0.8]
});

// Smooth lines
const smoothedLines = lines.map(line => ({
  ...line,
  paths: line.paths.map(path => chaikinSmooth(path, 3))
}));

// Generate and filter labels
const labels = removeCollidingLabels(
  generateLabelsForContours(smoothedLines),
  { width: 60, height: 20 }
);

// Render with labels
const svg = renderContoursSVG(
  { bands, lines: smoothedLines },
  { labels }
);
```

### Server-Side Rendering

```javascript
import { generateContoursAuto, renderContoursSVG } from 'd3-contour';
import { writeFile } from 'fs/promises';

const { bands, lines } = generateContoursAuto(grid, {
  size: [100, 100],
  count: 10
});

const svg = renderContoursSVG({ bands, lines }, {
  width: 1200,
  height: 900
});

await writeFile('output.svg', svg);
```

---

## Architecture

### Data Flow

```
Grid Data
  ↓
padGrid() (optional but recommended)
  ↓
generateContours() / generateContoursAuto()
  ↓
  { bands, lines }
  ↓
chaikinSmooth() / catmullRomSmooth() (on lines)
  ↓
generateLabelsForContours() → removeCollidingLabels()
  ↓
renderContoursSVG() / renderContoursCanvas()
  ↓
Output (SVG string or Canvas)
```

### Key Principles

1. **Isolines from Isobands**: Lines are extracted from band boundaries, not computed independently
2. **Post-Processing Smoothing**: Smooth the geometry, not the grid
3. **Unified Rendering**: Same geometry for both bands and lines
4. **Modular Design**: Use only what you need

---

## License

BSD-3-Clause (same as original d3-contour)
