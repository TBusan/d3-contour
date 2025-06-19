# Optimized d3-contour

An optimized version of the d3-contour library with improved performance, topology handling, edge clipping, and smoothing controls.

## Key Improvements

1. **Performance Optimization**
   - Added a spatial index for faster hole-polygon assignment
   - Replaced Array with Map for fragment tracking
   - Added bounds checking for early rejection in point-in-polygon tests
   - Optimized point-in-polygon algorithm

2. **Topology Correctness**
   - Implemented level-by-level contour processing to ensure topological consistency
   - Added post-processing to remove overlaps between different contour levels
   - Ensured unique attribution of areas to contour levels

3. **Edge Clipping**
   - Added proper edge clipping at grid boundaries
   - Preserved valid rings while discarding invalid ones after clipping

4. **Enhanced Smoothing Control**
   - Added a smoothFactor parameter to control the degree of smoothing
   - Preserved existing smooth/no-smooth option

## New API Methods

### Contours

```js
d3.contours()
    .size([width, height])
    .thresholds(thresholdValues)
    .smoothFactor(1.0)  // New: Control smoothing intensity (default: 1.0)
    .clipEdges(true)    // New: Enable/disable edge clipping (default: true)
```

### ContourDensity

```js
d3.contourDensity()
    .x(x)
    .y(y)
    .size([width, height])
    .cellSize(cellSize)
    .thresholds(thresholdValues)
    .bandwidth(bandwidth)
    .smoothFactor(1.0)  // New: Control smoothing intensity (default: 1.0)
    .clipEdges(true)    // New: Enable/disable edge clipping (default: true)
```

## Usage Example

```js
// Create contours with enhanced controls
const contours = d3.contours()
    .size([width, height])
    .thresholds(10)  // 10 contour levels
    .smoothFactor(0.8)  // Slightly reduced smoothing
    .clipEdges(true);  // Enable edge clipping

// Apply to data
const contourData = contours(values);

// Visualize
svg.selectAll("path")
    .data(contourData)
    .enter().append("path")
    .attr("d", d3.geoPath())
    .attr("fill", d => color(d.value));
```

## Technical Details

The optimizations focus on addressing key issues in the original implementation:

1. **Spatial Indexing**: Replaced O(n²) comparisons for hole assignments with a more efficient spatial index.

2. **Topological Consistency**: Implemented a level-by-level processing approach to ensure proper nesting relationships between contour levels.

3. **Edge Cases**: Added proper handling of grid edge cases with the new clipping functionality.

4. **Memory Efficiency**: Used Map objects for better performance with large datasets.

5. **Numerical Stability**: Improved collinearity tests for point-on-segment checks with a small epsilon value.

## Integration

This optimized library is a drop-in replacement for the original d3-contour library, with full compatibility with existing code while providing new optional features through additional parameters.
