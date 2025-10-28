# Enhanced D3-Contour

An enhanced version of d3-contour with significant improvements to address the limitations of the original implementation.

## Key Improvements

### 1. Saddle Point Disambiguation

The original d3-contour implementation had ambiguity in how to handle saddle points (cases 5 and 10 in the marching squares algorithm), which could lead to topological errors in the generated contours. This enhanced version:

- Adds proper saddle point disambiguation based on the average value at the center of each cell
- Ensures correct topology by choosing the appropriate contour configuration

### 2. Support for Multiple Output Modes

- **surfaces**: Generate contour polygons (MultiPolygon) - default mode
- **lines**: Generate contour lines (MultiLineString)
- **both**: Generate both representations simultaneously

### 3. Optimized Hole Assignment

The original implementation used an O(n²) algorithm for assigning holes to their containing polygons. This enhanced version:

- Uses a spatial index with bounding box optimization for faster containment tests
- Sorts candidates by area to ensure holes are assigned to the smallest containing polygon
- Significantly improves performance for complex contours with many holes

### 4. Enhanced Smoothing

- Configurable smoothing intensity instead of just on/off
- Improved algorithm that maintains the integrity of the contours while removing grid artifacts
- Better handling of boundary cases

### 5. Boundary Handling

- Improved handling of contours at grid boundaries
- Ability to close open contours at boundaries
- More accurate boundary detection and path generation

### 6. Self-Intersection Detection

- Added detection of topological errors like self-intersections
- Warning for self-intersecting contours to help identify problematic areas

### 7. Enhanced Kernel Density Estimation

- More accurate Gaussian kernel for density estimation
- Configurable kernel size and density factor
- Better handling of normalization
- Support for null value masking

### 8. GeoJSON Export

Added utility for converting contour output to GeoJSON format for easy visualization in mapping libraries.

## Usage Example

```javascript
import { contours, contourDensity, toGeoJSON } from "./src-ztContour2/index.js";

// Create basic contours for a grid of values
const data = [...]; // 1D array of values
const size = [width, height]; // Grid dimensions

// Create contour generator
const contour = contours()
  .size(size)
  .thresholds(10) // Either number of thresholds or array of values
  .mode("surfaces")
  .smooth(1.0);   // Smoothing level (0 = no smoothing)

// Generate contours
const result = contour(data);

// Convert to GeoJSON
const geoJson = toGeoJSON(result);

// For point density estimation
const density = contourDensity()
  .x(d => d.x)
  .y(d => d.y)
  .weight(d => d.value)
  .size([width, height])
  .bandwidth(20)
  .kernelSize(5)
  .kernelDensityFactor(1.5)
  .mode("surfaces");

// Generate density contours from point data
const densityContours = density(pointData);
```

## API Reference

### contours()

Creates a new contour generator with the default settings.

#### contours.size([width, height])
Sets the size of the grid.

#### contours.thresholds([count | values])
Sets the threshold values or count.

#### contours.smooth([level])
Sets the smoothing level (0 = off, higher values = more smoothing).

#### contours.mode(mode)
Sets the output mode ("surfaces", "lines", or "both").

#### contours.nullMask(mask)
Sets a mask array for null values.

### contourDensity()

Creates a new density contour generator.

#### density.x([accessor])
Sets the x-coordinate accessor.

#### density.y([accessor])
Sets the y-coordinate accessor.

#### density.weight([accessor])
Sets the point weight accessor.

#### density.size([width, height])
Sets the size of the output.

#### density.bandwidth([radius])
Sets the kernel bandwidth (blur radius).

#### density.cellSize([size])
Sets the cell size of the underlying grid.

#### density.kernelSize([size])
Sets the kernel size for density estimation.

#### density.kernelDensityFactor([factor])
Sets the density factor for the kernel.

#### density.thresholds([count | values])
Sets the threshold values or count.

#### density.mode(mode)
Sets the output mode ("surfaces", "lines", or "both").

#### density.nullValuesMask([boolean])
Enables/disables treating cells with zero density as null values.

### toGeoJSON(contours, [properties])

Converts contour output to GeoJSON Feature or FeatureCollection. 