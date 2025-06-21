# Custom Contour Library

This is a custom implementation of a contour and contour band generator based on the d3-contour library, with added support for:

1. Null values in the input data
2. Improved saddle point resolution for better contour topology
3. Contour bands (filled areas between contours)
4. Enhanced smoothing with boundary condition handling

## Usage

### Basic Contours Example

```js
import { contours } from './custom-src/index.js';

// Example data with null values
const data = [
  [1, 2, 3, 4, 5],
  [1, null, 3, 4, 5],
  [1, 2, 3, null, 5],
  [1, 2, 3, 4, 5],
  [0, 1, 2, 3, 4]
];

// Create contours
const contourGenerator = contours()
  .size([5, 5])       // Specify data dimensions
  .thresholds([1, 2, 3, 4]); // Specify threshold values

// Generate contours
const contourData = contourGenerator(data.flat());

// Each contour contains:
// - type: "MultiPolygon"
// - value: the threshold value
// - coordinates: GeoJSON MultiPolygon coordinates
```

### Contour Bands Example

```js
import { contourBands } from './custom-src/index.js';

// Create contour bands
const bandGenerator = contourBands()
  .size([5, 5])       // Specify data dimensions
  .thresholds([1, 2, 3, 4, 5]); // Specify threshold values

// Generate bands
const bandData = bandGenerator(data.flat());

// Each band contains:
// - type: "MultiPolygon"
// - lowerValue: lower threshold value
// - upperValue: upper threshold value
// - coordinates: GeoJSON MultiPolygon coordinates
```

## Key Features

### Null Value Handling

The library properly handles null, undefined, and NaN values in the input data by:
- Excluding them from threshold calculations
- Treating them as below threshold in the marching squares algorithm
- Handling them gracefully during smoothing operations

### Saddle Point Resolution

Saddle points can create ambiguity in the contour lines (cases 5 and 10 in marching squares).
This implementation resolves the ambiguity by looking at the average of the cell values to determine the correct configuration.

### Smooth Contours

The library includes enhanced smoothing that:
- Properly handles grid boundaries
- Gracefully manages null values
- Creates more visually appealing contours

### Contour Bands

The contour bands feature generates filled regions between threshold levels, useful for:
- Choropleth maps
- Heat maps
- Elevation visualization 