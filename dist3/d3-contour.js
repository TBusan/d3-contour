// https://d3js.org/d3-contour/ v4.0.2 Copyright 2012-2023 Mike Bostock
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports, require('d3-array')) :
typeof define === 'function' && define.amd ? define(['exports', 'd3-array'], factory) :
(global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.d3 = global.d3 || {}, global.d3));
})(this, (function (exports, d3Array) { 'use strict';

/**
 * Geometric utility functions for contour generation
 */

// Calculate the area of a ring (polygon)
// Positive for counter-clockwise (exterior rings), negative for clockwise (hole rings)
function area(ring) {
    let i = 0, n = ring.length;
    let area = ring[n - 1][1] * ring[0][0] - ring[n - 1][0] * ring[0][1];
    while (++i < n) area += ring[i - 1][1] * ring[i][0] - ring[i - 1][0] * ring[i][1];
    return area;
}

// Check if a ring contains a point
// Returns -1 for contained points, 0 for boundary points, 1 for exterior points
function contains(ring, hole) {
    let i = -1, n = hole.length, c;
    while (++i < n) if (c = ringContains(ring, hole[i])) return c;
    return 0;
}

// Helper function for contains() - tests if a single point is inside a ring
function ringContains(ring, point) {
    const [x, y] = point;
    let contains = -1;
    for (let i = 0, n = ring.length, j = n - 1; i < n; j = i++) {
        const [xi, yi] = ring[i];
        const [xj, yj] = ring[j];
        if (segmentContains(ring[i], ring[j], point)) return 0;
        if (((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
            contains = -contains;
        }
    }
    return contains;
}

// Check if a point is on a line segment
function segmentContains(a, b, c) {
    const i = +(a[0] === b[0]);
    return collinear(a, b, c) && within(a[i], c[i], b[i]);
}

// Check if three points are collinear (on the same line)
function collinear(a, b, c) {
    return (b[0] - a[0]) * (c[1] - a[1]) === (c[0] - a[0]) * (b[1] - a[1]);
}

// Check if a value is between two other values
function within(p, q, r) {
    return p <= q && q <= r || r <= q && q <= p;
}

// Check if two line segments intersect
function segmentIntersect(a, b, c, d) {
    // Calculate the direction vectors
    const ab = [b[0] - a[0], b[1] - a[1]];
    const cd = [d[0] - c[0], d[1] - c[1]];
    
    // Calculate the determinant
    const det = ab[0] * cd[1] - ab[1] * cd[0];
    
    // If determinant is zero, lines are parallel
    if (Math.abs(det) < 1e-10) return false;
    
    // Calculate the parameters t and s
    const ac = [c[0] - a[0], c[1] - a[1]];
    const t = (ac[0] * cd[1] - ac[1] * cd[0]) / det;
    const s = (ac[0] * ab[1] - ac[1] * ab[0]) / det;
    
    // Intersection occurs if t and s are both in [0,1]
    return t >= 0 && t <= 1 && s >= 0 && s <= 1;
}

/**
 * Optimized spatial index for polygon and hole matching
 * Addresses the O(n²) performance issue in the original implementation
 */


class SpatialIndex {
    constructor() {
        this.polygons = [];
        this.bounds = [];
    }
    
    // Add a polygon to the index with its bounding box
    addPolygon(polygon, index) {
        this.polygons.push({ polygon, index });
        this.bounds.push(this.computeBounds(polygon[0])); // Outer ring
    }
    
    // Calculate the bounding box of a ring
    computeBounds(ring) {
        let minX = Infinity, minY = Infinity;
        let maxX = -Infinity, maxY = -Infinity;
        
        for (const [x, y] of ring) {
            if (x < minX) minX = x;
            if (x > maxX) maxX = x;
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
        }
        
        return { minX, minY, maxX, maxY };
    }
    
    // Find the smallest polygon that contains the hole
    findContainingPolygon(hole) {
        const holeBounds = this.computeBounds(hole);
        const candidates = [];
        
        // First pass: bounding box test (much faster than full containment test)
        for (let i = 0; i < this.bounds.length; i++) {
            const bounds = this.bounds[i];
            if (holeBounds.minX >= bounds.minX && holeBounds.maxX <= bounds.maxX &&
                holeBounds.minY >= bounds.minY && holeBounds.maxY <= bounds.maxY) {
                candidates.push(i);
            }
        }
        
        // If we have no candidates after the bounding box test, return early
        if (candidates.length === 0) {
            return -1;
        }
        
        // For multiple candidates, we need to find the smallest containing polygon
        // This handles nested polygons correctly by preferring the innermost container
        let bestIndex = -1;
        let smallestArea = Infinity;
        
        for (const idx of candidates) {
            const polygon = this.polygons[idx].polygon;
            const ring = polygon[0];
            
            // Check if the hole is contained within this polygon
            if (contains(ring, hole) !== -1) {
                // Calculate approximate area to determine the smallest container
                const bounds = this.bounds[idx];
                const area = (bounds.maxX - bounds.minX) * (bounds.maxY - bounds.minY);
                
                if (area < smallestArea) {
                    smallestArea = area;
                    bestIndex = this.polygons[idx].index;
                }
            }
        }
        
        return bestIndex;
    }
}

/**
 * Helper functions for the enhanced contour library
 */

// No-operation function
function noop() {}

// When computing the extent, ignore infinite values (as well as invalid ones).
function finite(x) {
    return isFinite(x) ? x : NaN;
}

// Is the (possibly invalid) x greater than or equal to the (known valid) value?
// Treat any invalid value as below negative infinity.
function above(x, value) {
    return x != null && isFinite(x) && x >= value;
}

// During smoothing, treat any invalid value as negative infinity.
function valid(v) {
    return v == null || isNaN(v = +v) ? -Infinity : v;
}

// Check if two points are equal (within a small epsilon)
function pointsEqual(p1, p2) {
    const epsilon = 1e-6;
    return Math.abs(p1[0] - p2[0]) < epsilon && Math.abs(p1[1] - p2[1]) < epsilon;
}

/**
 * Topology utility functions for contour generation
 * Handles self-intersections and other topological errors
 */


// Check if a ring has self-intersections
function hasSelfIntersections(ring) {
    if (ring.length < 4) return false;
    
    // Check each pair of non-adjacent line segments for intersection
    for (let i = 0; i < ring.length - 1; i++) {
        const a = ring[i];
        const b = ring[i + 1];
        
        for (let j = i + 2; j < ring.length - 1; j++) {
            // Skip adjacent segments
            if (j === i - 1 || j === i || j === i + 1) continue;
            
            const c = ring[j];
            const d = ring[j + 1];
            
            // Skip if the segments share an endpoint
            if (pointsEqual(a, c) || pointsEqual(a, d) || pointsEqual(b, c) || pointsEqual(b, d)) continue;
            
            // Check for intersection
            if (segmentIntersect(a, b, c, d)) {
                return true;
            }
        }
    }
    
    return false;
}

// Remove self-intersections from a line
function removeSelfIntersections(line) {
    if (line.length < 4) return line;
    
    const result = [line[0]];
    let currentPoint = line[0];
    
    // Walk through the line, skipping points that would create self-intersections
    for (let i = 1; i < line.length; i++) {
        const nextPoint = line[i];
        let hasIntersection = false;
        
        // Check if adding this segment would create an intersection
        for (let j = 0; j < result.length - 1; j++) {
            if (segmentIntersect(
                currentPoint, nextPoint,
                result[j], result[j + 1]
            )) {
                hasIntersection = true;
                break;
            }
        }
        
        if (!hasIntersection) {
            result.push(nextPoint);
            currentPoint = nextPoint;
        }
    }
    
    return result;
}

/**
 * Enhanced smoothing functions for contour generation
 * Provides configurable smoothing with multiple strategies
 */


// Create a smoothing function with configurable parameters
function createSmoothing(smoothingLevel = 1.0) {
    if (smoothingLevel <= 0) return function() {}; // No smoothing
    
    return function smoothAdvanced(ring, values, value, dx, dy) {
        // Threshold for detecting grid-aligned points that need smoothing
        const SMOOTH_THRESHOLD = 0.03 * smoothingLevel; 
        
        ring.forEach(function(point) {
            const x = point[0];
            const y = point[1];
            const xt = Math.floor(x);
            const yt = Math.floor(y);
            
            if (xt >= 0 && xt < dx - 1 && yt >= 0 && yt < dy - 1) {
                // Get fractional position within grid cell
                const xFrac = x - xt;
                const yFrac = y - yt;
                
                // Detect and smooth grid-aligned points to avoid artifacts
                if (Math.abs(xFrac) < SMOOTH_THRESHOLD || Math.abs(xFrac - 1) < SMOOTH_THRESHOLD) {
                    // Smoothly adjust x coordinate away from grid lines
                    point[0] += (xFrac < 0.5 ? 1 : -1) * SMOOTH_THRESHOLD * smoothingLevel * 0.5;
                }
                
                if (Math.abs(yFrac) < SMOOTH_THRESHOLD || Math.abs(yFrac - 1) < SMOOTH_THRESHOLD) {
                    // Smoothly adjust y coordinate away from grid lines
                    point[1] += (yFrac < 0.5 ? 1 : -1) * SMOOTH_THRESHOLD * smoothingLevel * 0.5;
                }
                
                // Enhanced interpolation at cell edges
                if (Math.abs(xFrac - 0.5) < SMOOTH_THRESHOLD && x > 0 && x < dx) {
                    const v0 = valid(values[yt * dx + xt]);
                    const v1 = valid(values[yt * dx + xt + 1]);
                    point[0] = smoothInterpolate(x, v0, v1, value, smoothingLevel);
                }
                
                if (Math.abs(yFrac - 0.5) < SMOOTH_THRESHOLD && y > 0 && y < dy) {
                    const v0 = valid(values[yt * dx + xt]);
                    const v1 = valid(values[(yt + 1) * dx + xt]);
                    point[1] = smoothInterpolate(y, v0, v1, value, smoothingLevel);
                }
            }
        });
    };
}

// Enhanced interpolation function with smoothing factor
function smoothInterpolate(coord, v0, v1, value, smoothingFactor = 1.0) {
    // Get the base coordinate (integer part)
    const base = Math.floor(coord);
    
    // Handle edge cases
    if (v0 === v1) return base + 0.5;
    if (!isFinite(v0) || !isFinite(v1)) return coord;
    
    // Calculate interpolation factor
    const a = value - v0;
    const b = v1 - v0;
    let d = isFinite(a) && isFinite(b) && b !== 0 ? a / b : 0.5;
    
    // Apply smoothing factor - higher values make interpolation more aggressive
    if (smoothingFactor !== 1.0) {
        // Move interpolation factor toward 0.5 for smoother transitions
        d = d + (0.5 - d) * (1 - smoothingFactor);
    }
    
    // Keep within valid range
    if (d < 0) d = 0;
    if (d > 1) d = 1;
    
    // Apply simple linear interpolation
    return base + d;
}

/**
 * Enhanced contours implementation addressing limitations in the original d3-contour
 * 
 * Improvements:
 * 1. Proper saddle point disambiguation
 * 2. Support for both contour lines and surfaces
 * 3. Optimized hole assignment with spatial indexing
 * 4. Enhanced smoothing with configurable parameters
 * 5. Null value masking
 * 6. Boundary handling
 */


// Enhanced marching squares cases with proper orientation
// Each case is an array of line segments
// Each segment connects two points on the edges of a unit square
const marchingSquaresCases = [
    [],                                     // Case 0: No contour
    [[[0, 0.5], [0.5, 0]]],                // Case 1: Bottom-left corner
    [[[0.5, 0], [1, 0.5]]],                // Case 2: Bottom-right corner
    [[[0, 0.5], [1, 0.5]]],                // Case 3: Bottom edge
    [[[1, 0.5], [0.5, 1]]],                // Case 4: Top-right corner
    [[[0, 0.5], [0.5, 0]], [[1, 0.5], [0.5, 1]]], // Case 5: Saddle (bottom-left + top-right)
    [[[0.5, 0], [0.5, 1]]],                // Case 6: Right edge
    [[[0, 0.5], [0.5, 1], [0.5, 0]]],      // Case 7: Not bottom-right corner
    [[[0.5, 1], [0, 0.5]]],                // Case 8: Top-left corner
    [[[0.5, 0], [0.5, 1]]],                // Case 9: Left edge
    [[[0.5, 1], [0, 0.5]], [[0.5, 0], [1, 0.5]]], // Case 10: Saddle (top-left + bottom-right)
    [[[0.5, 1], [1, 0.5], [0.5, 0]]],      // Case 11: Not bottom-left corner
    [[[0.5, 1], [0, 0.5], [1, 0.5]]],      // Case 12: Top edge
    [[[0.5, 0], [0, 0.5], [0.5, 1]]],      // Case 13: Not top-right corner
    [[[0.5, 0], [1, 0.5], [0.5, 1]]],      // Case 14: Not top-left corner
    []                                      // Case 15: No contour
];

// Enhanced contour generator
function contours() {
    let dx = 1;
    let dy = 1;
    let threshold = d3Array.thresholdSturges;
    let smoothingLevel = 1.0;
    let smooth = createSmoothing(smoothingLevel);
    let mode = "surfaces"; // "surfaces", "lines", or "both"
    let nullMask = null;

    function contours(values) {
        const tz = computeThresholds(values);
        
        if (mode === "lines") {
            return tz.map(value => contourLines(values, value));
        } else if (mode === "surfaces") {
            return tz.map(value => contourSurfaces(values, value));
        } else { // both
            return tz.map(value => ({
                lines: contourLines(values, value),
                surfaces: contourSurfaces(values, value),
                value: value
            }));
        }
    }
    
    function computeThresholds(values) {
        let tz = threshold(values);
        
        // Convert number of thresholds into uniform thresholds
        if (!Array.isArray(tz)) {
            const e = d3Array.extent(values, finite);
            if (isNaN(e[0]) || isNaN(e[1])) return [];
            tz = d3Array.ticks(e[0], e[1], tz);
        } else {
            tz = tz.slice().sort((a, b) => a - b);
        }
        
        return tz;
    }
    
    // Generate contour lines (open paths)
    function contourLines(values, value) {
        const v = value == null ? NaN : +value;
        if (isNaN(v)) throw new Error(`invalid value: ${value}`);
        
        const lines = [];
        
        // Use marching squares to generate line segments
        marchingSquares(values, v, (x, y, caseIndex, corners) => {
            const segments = getContourSegments(caseIndex, corners, v);
            
            // Convert segments to absolute coordinates
            segments.forEach(segment => {
                if (segment.length < 2) return; // Skip invalid segments
                
                const line = segment.map(point => [
                    point[0] + x,
                    point[1] + y
                ]);
                
                // Apply smoothing if enabled
                if (smooth) {
                    smooth(line, values, v, dx, dy);
                }
                
                lines.push(line);
            });
        });
        
        // Connect line segments into continuous lines
        let connectedLines = connectSegments(lines);
        
        // Filter out very short lines and grid artifacts
        connectedLines = connectedLines.filter(line => {
            // Remove short lines (likely artifacts)
            if (line.length < 4) return false;
            
            // Calculate the path length
            let pathLength = 0;
            for (let i = 1; i < line.length; i++) {
                pathLength += distance(line[i-1], line[i]);
            }
            
            // Filter out lines that are too short - REDUCED THRESHOLD FROM 2.0 to 0.5
            if (pathLength < 0.5) return false;
            
            // Filter out grid-aligned artifacts - RELAXED THRESHOLD FROM 0.5 to 0.8
            let verticalSegments = 0;
            let horizontalSegments = 0;
            
            for (let i = 1; i < line.length; i++) {
                const dx = Math.abs(line[i][0] - line[i-1][0]);
                const dy = Math.abs(line[i][1] - line[i-1][1]);
                
                if (dx < 1e-4) verticalSegments++;
                if (dy < 1e-4) horizontalSegments++;
            }
            
            // If more than 80% (instead of 50%) of segments are grid-aligned, it's likely an artifact
            const totalSegments = line.length - 1;
            if (verticalSegments > totalSegments * 0.8 || horizontalSegments > totalSegments * 0.8) {
                return false;
            }
            
            return true;
        });
        
        // Try to close open loops that are nearly closed
        connectedLines = connectedLines.map(line => {
            if (line.length >= 4) {
                const start = line[0];
                const end = line[line.length - 1];
                const dist = distance(start, end);
                
                // If the start and end are reasonably close, close the loop
                if (dist < 2.0 && dist > 1e-4) {
                    return [...line, line[0].slice()];
                }
            }
            return line;
        });
        
        // Apply additional smoothing to reduce sharp angles and self-intersections
        connectedLines = smoothLines(connectedLines);
        connectedLines = connectedLines.map(removeSelfIntersections);
        
        return {
            type: "MultiLineString",
            value: value,
            coordinates: connectedLines
        };
    }
    
    // Generate contour surfaces (polygons) with enhanced hole assignment
    function contourSurfaces(values, value) {
        const v = value == null ? NaN : +value;
        if (isNaN(v)) throw new Error(`invalid value: ${value}`);
        
        // First get the contour lines
        const lines = [];
        
        // Use marching squares to generate line segments
        marchingSquares(values, v, (x, y, caseIndex, corners) => {
            const segments = getContourSegments(caseIndex, corners, v);
            
            // Convert segments to absolute coordinates
            segments.forEach(segment => {
                const line = segment.map(point => [
                    point[0] + x,
                    point[1] + y
                ]);
                
                // Apply smoothing if enabled
                if (smooth) {
                    smooth(line, values, v, dx, dy);
                }
                
                lines.push(line);
            });
        });
        
        // Connect line segments into closed rings
        const rings = connectSegmentsToRings(lines);
        
        if (rings.length === 0) {
            return {
                type: "MultiPolygon",
                value: value,
                coordinates: []
            };
        }
        
        // Separate rings into polygons and holes based on winding order
        const polygons = [];
        const holes = [];
        
        rings.forEach(ring => {
            // Skip invalid rings (too small or self-intersecting)
            if (ring.length < 4 || hasSelfIntersections(ring)) return;
            
            // Calculate area to determine if it's a hole or exterior
            const ringArea = area(ring);
            
            if (ringArea > 0) {
                // Counter-clockwise rings are polygon exteriors
                polygons.push([ring]);
            } else if (ringArea < 0) {
                // Clockwise rings are holes
                holes.push(ring);
            }
        });
        
        // Use an optimized spatial index for hole assignment
        const spatialIndex = new SpatialIndex();
        polygons.forEach((polygon, i) => spatialIndex.addPolygon(polygon, i));
        
        // Assign holes to their containing polygons using the spatial index
        holes.forEach(hole => {
            const containingIndex = spatialIndex.findContainingPolygon(hole);
            if (containingIndex !== -1) {
                polygons[containingIndex].push(hole);
            }
        });
        
        // Filter out invalid polygons (e.g., with self-intersections)
        const validPolygons = polygons.filter(polygon => {
            // Check if the exterior ring is valid
            return polygon[0].length >= 4;
        });
        
        return {
            type: "MultiPolygon",
            value: value,
            coordinates: validPolygons
        };
    }
    
    // Core marching squares algorithm with saddle point disambiguation
    function marchingSquares(values, value, callback) {
        // Process each cell in the grid
        for (let y = 0; y < dy - 1; y++) {
            for (let x = 0; x < dx - 1; x++) {
                // Get the values at the four corners of the cell
                const corners = [
                    [getValue(values, x, y), getValue(values, x + 1, y)],
                    [getValue(values, x, y + 1), getValue(values, x + 1, y + 1)]
                ];
                
                // LESS STRICT: Handle NaN values by trying to interpolate
                let hasNaN = false;
                for (let i = 0; i < 2; i++) {
                    for (let j = 0; j < 2; j++) {
                        if (isNaN(corners[i][j])) {
                            hasNaN = true;
                            // Try to interpolate from neighbors
                            corners[i][j] = interpolateNaN(values, x + j, y + i, dx, dy);
                        }
                    }
                }
                
                // Skip cells with too many NaN values or if interpolation wasn't successful
                if ((hasNaN && (isNaN(corners[0][0]) || isNaN(corners[0][1]) || 
                    isNaN(corners[1][0]) || isNaN(corners[1][1]))) ||
                    (isNaN(corners[0][0]) && isNaN(corners[0][1]) && 
                    isNaN(corners[1][0]) && isNaN(corners[1][1]))) {
                    continue;
                }
                
                // Determine the case index (0-15) for this cell
                const caseIndex = getMarchingIndex(value, corners);
                
                // Skip cases with no contour lines
                if (caseIndex === 0 || caseIndex === 15) continue;
                
                // Handle saddle cases with proper disambiguation
                const finalCaseIndex = disambiguateSaddle(caseIndex, corners, value);
                
                // Call the callback with the cell info
                callback(x, y, finalCaseIndex, corners);
            }
        }
    }
    
    // Helper function to interpolate NaN values from neighbors
    function interpolateNaN(values, x, y, dx, dy) {
        const neighbors = [];
        
        // Check 8 neighboring cells
        for (let ny = Math.max(0, y - 1); ny <= Math.min(dy - 1, y + 1); ny++) {
            for (let nx = Math.max(0, x - 1); nx <= Math.min(dx - 1, x + 1); nx++) {
                if (nx === x && ny === y) continue; // Skip self
                
                const val = getValue(values, nx, ny);
                if (!isNaN(val)) {
                    neighbors.push(val);
                }
            }
        }
        
        // If we have neighbors, average them
        if (neighbors.length > 0) {
            return neighbors.reduce((a, b) => a + b, 0) / neighbors.length;
        }
        
        // Otherwise, return NaN
        return NaN;
    }
    
    function getValue(values, x, y) {
        if (x < 0 || x >= dx || y < 0 || y >= dy) return NaN;
        const index = y * dx + x;
        if (index < 0 || index >= values.length) return NaN;
        const value = values[index];
        // LESS STRICT: Only mask if nullMask is provided and the value is explicitly 0
        return nullMask && nullMask[index] === 0 ? NaN : value;
    }
    
    // Enhanced saddle disambiguation (inspired by plotly.js)
    function disambiguateSaddle(caseIndex, corners, value) {
        if (caseIndex !== 5 && caseIndex !== 10) return caseIndex;
        
        // Calculate average value at center of the cell
        const avg = (corners[0][0] + corners[0][1] + corners[1][0] + corners[1][1]) / 4;
        
        // Determine which configuration to use based on the center value
        if (caseIndex === 5) {
            // Two possible configurations for case 5
            return value > avg ? 5 : 10;
        } else { // caseIndex === 10
            // Two possible configurations for case 10
            return value > avg ? 10 : 5;
        }
    }
    
    function getMarchingIndex(value, corners) {
        // Calculate index based on which corners are above the threshold
        return (above(corners[0][0], value) ? 1 : 0) +
               (above(corners[0][1], value) ? 2 : 0) +
               (above(corners[1][1], value) ? 4 : 0) +
               (above(corners[1][0], value) ? 8 : 0);
    }
    
    // Get contour line segments for a given case
    function getContourSegments(caseIndex, corners, value) {
        // Get the basic segments for this case
        const segments = marchingSquaresCases[caseIndex];
        
        if (!segments || segments.length === 0) {
            return [];
        }
        
        // For each segment, interpolate the exact crossing points
        return segments.map(segment => {
            return segment.map(point => {
                const [x, y] = point;
                
                // If the point is already at a corner, no interpolation needed
                if ((x === 0 || x === 1) && (y === 0 || y === 1)) {
                    return point;
                }
                
                // Interpolate the exact crossing point
                if (x === 0.5) {
                    // Crossing on vertical edge
                    const y0 = y === 0 ? 0 : 1;
                    const v0 = corners[y0][0];
                    const v1 = corners[y0][1];
                    const t = linearInterpolation(v0, v1, value);
                    return [t, y];
                } else if (y === 0.5) {
                    // Crossing on horizontal edge
                    const x0 = x === 0 ? 0 : 1;
                    const v0 = corners[0][x0];
                    const v1 = corners[1][x0];
                    const t = linearInterpolation(v0, v1, value);
                    return [x, t];
                }
                
                // Should never reach here
                return point;
            });
        });
    }
    
    // Linear interpolation helper
    function linearInterpolation(v0, v1, value) {
        // Handle edge cases
        if (v0 === v1) return 0.5;
        if (!isFinite(v0) || !isFinite(v1)) return 0.5;
        
        // Calculate interpolation factor
        const t = (value - v0) / (v1 - v0);
        
        // Ensure result is within [0,1] range and avoid numerical issues
        if (!isFinite(t)) return 0.5;
        if (t < 0) return 0;
        if (t > 1) return 1;
        return t;
    }
    
    // Connect line segments into continuous lines
    function connectSegments(segments) {
        if (segments.length === 0) return [];
        
        // Create a deep copy of segments to avoid modifying the original
        const segmentsCopy = segments.map(segment => segment.slice());
        
        // Sort segments by their starting point's x-coordinate for more consistent joining
        segmentsCopy.sort((a, b) => a[0][0] - b[0][0]);
        
        const result = [];
        const used = new Set();
        
        // Start with any unused segment
        for (let i = 0; i < segmentsCopy.length; i++) {
            if (used.has(i)) continue;
            
            used.add(i);
            let line = segmentsCopy[i].slice();
            
            // Try to extend the line in both directions
            let extended;
            do {
                extended = false;
                
                // Try to extend at the end
                const end = line[line.length - 1];
                let bestMatch = -1;
                // INCREASED threshold for matching points (from 1e-4 to 1e-3)
                let bestMatchDistance = 1e-3; 
                
                for (let j = 0; j < segmentsCopy.length; j++) {
                    if (used.has(j)) continue;
                    
                    const segment = segmentsCopy[j];
                    const start = segment[0];
                    const last = segment[segment.length - 1];
                    
                    // Calculate distances
                    const distToStart = distance(end, start);
                    const distToEnd = distance(end, last);
                    
                    // Find the best match (closest point)
                    if (distToStart < bestMatchDistance) {
                        bestMatch = j;
                        bestMatchDistance = distToStart;
                        // Flag for connecting to start of segment
                        segment._connectToStart = true;
                    } else if (distToEnd < bestMatchDistance) {
                        bestMatch = j;
                        bestMatchDistance = distToEnd;
                        // Flag for connecting to end of segment
                        segment._connectToStart = false;
                    }
                }
                
                // If we found a good match, extend the line
                if (bestMatch !== -1) {
                    const segment = segmentsCopy[bestMatch];
                    if (segment._connectToStart) {
                        // Add all points except the first (which is close to our end)
                        line.push(...segment.slice(1));
                    } else {
                        // Add all points except the last (which is close to our end), in reverse order
                        line.push(...segment.slice(0, -1).reverse());
                    }
                    used.add(bestMatch);
                    extended = true;
                }
                
                if (!extended) {
                    // Try to extend at the start
                    const start = line[0];
                    bestMatch = -1;
                    bestMatchDistance = 1e-3;  // INCREASED threshold for matching points
                    
                    for (let j = 0; j < segmentsCopy.length; j++) {
                        if (used.has(j)) continue;
                        
                        const segment = segmentsCopy[j];
                        const segStart = segment[0];
                        const segEnd = segment[segment.length - 1];
                        
                        // Calculate distances
                        const distToStart = distance(start, segStart);
                        const distToEnd = distance(start, segEnd);
                        
                        // Find the best match (closest point)
                        if (distToStart < bestMatchDistance) {
                            bestMatch = j;
                            bestMatchDistance = distToStart;
                            segment._connectToStart = true;
                        } else if (distToEnd < bestMatchDistance) {
                            bestMatch = j;
                            bestMatchDistance = distToEnd;
                            segment._connectToStart = false;
                        }
                    }
                    
                    // If we found a good match, extend the line
                    if (bestMatch !== -1) {
                        const segment = segmentsCopy[bestMatch];
                        if (segment._connectToStart) {
                            // Add all points except the first (which is close to our start), in reverse order
                            line.unshift(...segment.slice(1).reverse());
                        } else {
                            // Add all points except the last (which is close to our start)
                            line.unshift(...segment.slice(0, -1));
                        }
                        used.add(bestMatch);
                        extended = true;
                    }
                }
            } while (extended);
            
            // Check if the line forms a closed loop
            if (line.length > 3 && distance(line[0], line[line.length - 1]) < 1e-3) {  // INCREASED threshold
                // Make sure the first and last points are exactly the same for a proper closed loop
                line[line.length - 1] = line[0].slice();
            }
            
            // Add the completed line to the result
            result.push(line);
        }
        
        return result;
    }
    
    // Connect line segments into closed rings
    function connectSegmentsToRings(segments) {
        if (segments.length === 0) return [];
        
        // First, connect segments into lines
        const lines = connectSegments(segments);
        const rings = [];
        
        // Process each line
        lines.forEach(line => {
            if (line.length < 3) return; // Skip too short lines
            
            const start = line[0];
            const end = line[line.length - 1];
            
            // Check if the line is already a ring (closed loop)
            if (distance(start, end) < 1e-4) {
                // Ensure first and last points are exactly the same
                const closedRing = line.slice(0, -1);
                closedRing.push(closedRing[0].slice());
                rings.push(closedRing);
                return;
            }
            
            // Handle boundary cases - try to close the ring along the boundary
            // This is for contours that intersect the grid boundary
            if (isBoundaryPoint(start, dx, dy) && isBoundaryPoint(end, dx, dy)) {
                // Create a closed ring by adding boundary points
                const closedRing = createBoundaryClosure(line, dx, dy);
                if (closedRing && closedRing.length >= 3) {
                    rings.push(closedRing);
                    return;
                }
            }
            
            // If we can't form a proper ring, don't include this line
            // This helps prevent topological errors in the surfaces
        });
        
        return rings;
    }
    
    // Helper function to calculate distance between two points
    function distance(p1, p2) {
        const dx = p1[0] - p2[0];
        const dy = p1[1] - p2[1];
        return Math.sqrt(dx * dx + dy * dy);
    }
    
    // Check if a point is on the boundary of the grid
    function isBoundaryPoint(point, dx, dy) {
        const [x, y] = point;
        const epsilon = 1e-4;
        return x < epsilon || x > dx - epsilon || y < epsilon || y > dy - epsilon;
    }
    
    // Create a closure for an open line along the boundary
    function createBoundaryClosure(line, dx, dy) {
        const start = line[0];
        const end = line[line.length - 1];
        
        // Don't try to close if points are too far apart
        if (distance(start, end) > Math.max(dx, dy) / 2) {
            return null;
        }
        
        const result = line.slice();
        const boundaryPoints = [];
        
        // Determine which boundaries the points are on
        const startOnLeft = start[0] < 1e-4;
        const startOnRight = start[0] > dx - 1e-4;
        const startOnBottom = start[1] < 1e-4;
        const startOnTop = start[1] > dy - 1e-4;
        
        const endOnLeft = end[0] < 1e-4;
        const endOnRight = end[0] > dx - 1e-4;
        const endOnBottom = end[1] < 1e-4;
        const endOnTop = end[1] > dy - 1e-4;
        
        // Simple case: both points on the same boundary
        if ((startOnLeft && endOnLeft) || 
            (startOnRight && endOnRight) || 
            (startOnTop && endOnTop) || 
            (startOnBottom && endOnBottom)) {
            // Just connect directly
            result.push(result[0].slice());
            return result;
        }
        
        // Corner cases
        if ((startOnLeft && endOnBottom) || (startOnBottom && endOnLeft)) {
            boundaryPoints.push([0, 0]); // Bottom-left corner
        } else if ((startOnRight && endOnBottom) || (startOnBottom && endOnRight)) {
            boundaryPoints.push([dx, 0]); // Bottom-right corner
        } else if ((startOnLeft && endOnTop) || (startOnTop && endOnLeft)) {
            boundaryPoints.push([0, dy]); // Top-left corner
        } else if ((startOnRight && endOnTop) || (startOnTop && endOnRight)) {
            boundaryPoints.push([dx, dy]); // Top-right corner
        } else {
            // More complex boundary traversal
            if (startOnLeft && endOnRight) {
                // Go through bottom
                boundaryPoints.push([0, 0], [dx, 0]);
            } else if (startOnRight && endOnLeft) {
                // Go through top
                boundaryPoints.push([dx, dy], [0, dy]);
            } else if (startOnBottom && endOnTop) {
                // Go through right
                boundaryPoints.push([dx, 0], [dx, dy]);
            } else if (startOnTop && endOnBottom) {
                // Go through left
                boundaryPoints.push([0, dy], [0, 0]);
            }
        }
        
        // Add boundary points and close the loop
        if (boundaryPoints.length > 0) {
            result.push(...boundaryPoints);
            result.push(result[0].slice());
            return result;
        }
        
        return null;
    }
    
    // New function to smooth lines to remove sharp angles
    function smoothLines(lines) {
        return lines.map(line => {
            if (line.length < 4) return line;
            
            const smoothed = [line[0]];
            
            // Use a simple moving average for smoothing
            for (let i = 1; i < line.length - 1; i++) {
                const prev = line[i - 1];
                const curr = line[i];
                const next = line[i + 1];
                
                // Skip points that create sharp angles
                const v1 = [curr[0] - prev[0], curr[1] - prev[1]];
                const v2 = [next[0] - curr[0], next[1] - curr[1]];
                
                const len1 = Math.sqrt(v1[0]*v1[0] + v1[1]*v1[1]);
                const len2 = Math.sqrt(v2[0]*v2[0] + v2[1]*v2[1]);
                
                // Skip if either segment is very short
                if (len1 < 1e-4 || len2 < 1e-4) continue;
                
                // Normalize vectors
                v1[0] /= len1; v1[1] /= len1;
                v2[0] /= len2; v2[1] /= len2;
                
                // Calculate dot product
                const dotProduct = v1[0]*v2[0] + v1[1]*v2[1];
                
                // If angle is too sharp (dot product too negative), smooth the point
                if (dotProduct < -0.7) {
                    // Create a smoothed point
                    const smoothedPoint = [
                        (prev[0] + curr[0] + next[0]) / 3,
                        (prev[1] + curr[1] + next[1]) / 3
                    ];
                    smoothed.push(smoothedPoint);
                } else {
                    smoothed.push(curr);
                }
            }
            
            // Add the last point
            if (line.length > 1) {
                smoothed.push(line[line.length - 1]);
            }
            
            return smoothed;
        });
    }
    
    // API Methods
    contours.contour = function(values, value) {
        if (mode === "lines") {
            return contourLines(values, value);
        } else if (mode === "both") {
            return {
                lines: contourLines(values, value),
                surfaces: contourSurfaces(values, value),
                value: value
            };
        }
        return contourSurfaces(values, value);
    };
    
    contours.size = function(_) {
        if (!arguments.length) return [dx, dy];
        const _0 = Math.floor(_[0]), _1 = Math.floor(_[1]);
        if (!(_0 >= 0 && _1 >= 0)) throw new Error("invalid size");
        return dx = _0, dy = _1, contours;
    };
    
    contours.thresholds = function(_) {
        return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? () => _.slice() : () => _, contours) : threshold;
    };
    
    contours.smooth = function(_) {
        if (!arguments.length) return smoothingLevel;
        if (typeof _ === "number") {
            smoothingLevel = _;
            smooth = createSmoothing(_);
        } else {
            smoothingLevel = _ ? 1.0 : 0;
            smooth = _ ? createSmoothing(1.0) : noop;
        }
        return contours;
    };
    
    contours.mode = function(_) {
        if (!arguments.length) return mode;
        if (!["lines", "surfaces", "both"].includes(_)) throw new Error("invalid mode");
        mode = _;
        return contours;
    };
    
    contours.nullMask = function(_) {
        if (!arguments.length) return nullMask;
        nullMask = _;
        return contours;
    };
    
    return contours;
}

/**
 * GeoJSON export utilities for contour generation
 * Provides tools to convert contour data to standard GeoJSON formats
 */

// Convert contour output to GeoJSON feature or feature collection
// Works with lines, surfaces or both
function toGeoJSON(contourResult, properties = {}) {
    if (Array.isArray(contourResult)) {
        return {
            type: "FeatureCollection",
            features: contourResult.map((item, i) => toGeoJSON(item, { ...properties, index: i }))
        };
    }
    
    if (contourResult.lines && contourResult.surfaces) {
        // Both mode
        return {
            type: "FeatureCollection",
            features: [
                {
                    type: "Feature",
                    properties: { ...properties, type: "lines", value: contourResult.value },
                    geometry: contourResult.lines
                },
                {
                    type: "Feature",
                    properties: { ...properties, type: "surfaces", value: contourResult.value },
                    geometry: contourResult.surfaces
                }
            ]
        };
    }
    
    return {
        type: "Feature",
        properties: { ...properties, value: contourResult.value },
        geometry: contourResult
    };
}

// Apply a projection function to contour geometry
function transformContours(contourResult, projectionFn) {
    if (Array.isArray(contourResult)) {
        return contourResult.map(item => transformContours(item, projectionFn));
    }
    
    if (contourResult.lines && contourResult.surfaces) {
        return {
            lines: transformGeometry(contourResult.lines, projectionFn),
            surfaces: transformGeometry(contourResult.surfaces, projectionFn),
            value: contourResult.value
        };
    }
    
    return transformGeometry(contourResult, projectionFn);
}

// Transform a GeoJSON geometry with a projection function
function transformGeometry(geometry, projectionFn) {
    if (!geometry || !geometry.coordinates) return geometry;
    
    const transformed = {
        ...geometry,
        coordinates: transformCoordinates(geometry.coordinates, projectionFn)
    };
    
    return transformed;
}

// Transform coordinates array recursively
function transformCoordinates(coordinates, projectionFn) {
    if (!Array.isArray(coordinates)) return coordinates;
    
    if (Array.isArray(coordinates[0])) {
        if (typeof coordinates[0][0] === 'number' && coordinates[0].length >= 2) {
            // This is a point, apply projection
            return projectionFn(coordinates);
        }
        // Recurse into nested arrays
        return coordinates.map(coord => transformCoordinates(coord, projectionFn));
    }
    
    return coordinates;
}

/**
 * Enhanced contour density generator
 * Improves the original d3-contour density with better handling of smoothing and nulls
 */


// Default accessors
function defaultX(d) {
    return d[0];
}

function defaultY(d) {
    return d[1];
}

function defaultWeight(d) {
    return 1;
}

// Enhanced density contour generator
function density() {
    let x = defaultX,
        y = defaultY,
        weight = defaultWeight,
        dx = 960,
        dy = 500,
        r = 20, // blur radius
        k = 2, // log2(grid cell size)
        o = r * 3, // grid offset, to pad for blur
        n = (dx + o * 2) >> k, // grid width
        m = (dy + o * 2) >> k, // grid height
        threshold = function() { return 20; },
        mode = "surfaces",
        bandwidthAdjust = 1.0,
        nullValuesMask = false;
    
    // Generate the grid for density estimation
    function grid(data) {
        // Use typed array for better performance
        var values = new Float32Array(n * m),
            pow2k = Math.pow(2, -k),
            i = -1;
        
        // Initialize grid with zeros
        for (let j = 0; j < values.length; j++) {
            values[j] = 0;
        }
        
        // Accumulate point weights to grid cells
        for (const d of data) {
            var xi = (x(d, ++i, data) + o) * pow2k,
                yi = (y(d, i, data) + o) * pow2k,
                wi = +weight(d, i, data);
            
            if (wi && xi >= 0 && xi < n && yi >= 0 && yi < m) {
                var x0 = Math.floor(xi),
                    y0 = Math.floor(yi),
                    xt = xi - x0 - 0.5,
                    yt = yi - y0 - 0.5;
                
                // Bilinear interpolation to distribute weight to surrounding cells
                values[x0 + y0 * n] += (1 - xt) * (1 - yt) * wi;
                values[x0 + 1 + y0 * n] += xt * (1 - yt) * wi;
                values[x0 + 1 + (y0 + 1) * n] += xt * yt * wi;
                values[x0 + (y0 + 1) * n] += (1 - xt) * yt * wi;
            }
        }
        
        // Apply blur with adjustable radius
        d3Array.blur2({data: values, width: n, height: m}, r * pow2k * bandwidthAdjust);
        return values;
    }
    
    // Generate contours from density grid
    function density(data) {
        var values = grid(data),
            tz = threshold(values),
            pow4k = Math.pow(2, 2 * k),
            nulls = nullValuesMask ? createNullMask(values) : null;
        
        // Convert number of thresholds into uniform thresholds
        if (!Array.isArray(tz)) {
            // Find non-zero maximum for better thresholding
            const maxVal = d3Array.max(values) / pow4k;
            if (maxVal <= 0) return []; // No data
            
            // Generate thresholds
            tz = d3Array.ticks(Number.MIN_VALUE, maxVal, tz);
        }
        
        // Generate contours with the specified mode
        const contourGen = contours()
            .size([n, m])
            .thresholds(tz.map(function(threshold) { return threshold * pow4k; }))
            .mode(mode);
        
        // Apply null mask if enabled
        if (nulls) {
            contourGen.nullMask(nulls);
        }
        
        // Generate and transform contours
        return contourGen(values).map((c, i) => {
            if (mode === "both") {
                c.lines = transform(c.lines);
                c.surfaces = transform(c.surfaces);
                c.value = +tz[i];
                return c;
            } else {
                const transformed = transform(c);
                transformed.value = +tz[i];
                return transformed;
            }
        });
    }
    
    // Create a special contour function for individual values
    density.contours = function(data) {
        var values = grid(data),
            contourGen = contours().size([n, m]).mode(mode),
            pow4k = Math.pow(2, 2 * k),
            nulls = nullValuesMask ? createNullMask(values) : null;
        
        // Apply null mask if enabled
        if (nulls) {
            contourGen.nullMask(nulls);
        }
        
        // Function to generate contours for a specific value
        const contour = value => {
            value = +value;
            
            if (mode === "both") {
                const result = contourGen.contour(values, value * pow4k);
                result.lines = transform(result.lines);
                result.surfaces = transform(result.surfaces);
                result.value = value; // preserve exact value
                return result;
            } else {
                const c = transform(contourGen.contour(values, value * pow4k));
                c.value = value; // preserve exact value
                return c;
            }
        };
        
        // Add max value property for convenience
        Object.defineProperty(contour, "max", {
            get: () => d3Array.max(values) / pow4k
        });
        
        return contour;
    };
    
    // Create mask for null values
    function createNullMask(values) {
        const mask = new Uint8Array(values.length);
        for (let i = 0; i < values.length; i++) {
            // Mark cells with zero density as null
            mask[i] = values[i] > 0 ? 1 : 0;
        }
        return mask;
    }
    
    // Transform coordinates from grid to original space
    function transform(geometry) {
        return transformContours(geometry, transformPoint);
    }
    
    // Transform a single point from grid to original space
    function transformPoint(point) {
        // Clone the point to avoid modifying the original
        const result = [0, 0];
        result[0] = point[0] * Math.pow(2, k) - o;
        result[1] = point[1] * Math.pow(2, k) - o;
        return result;
    }
    
    // API methods
    density.x = function(_) {
        return arguments.length ? (x = typeof _ === "function" ? _ : () => +_, density) : x;
    };
    
    density.y = function(_) {
        return arguments.length ? (y = typeof _ === "function" ? _ : () => +_, density) : y;
    };
    
    density.weight = function(_) {
        return arguments.length ? (weight = typeof _ === "function" ? _ : () => +_, density) : weight;
    };
    
    density.size = function(_) {
        if (!arguments.length) return [dx, dy];
        var _0 = +_[0], _1 = +_[1];
        if (!(_0 >= 0 && _1 >= 0)) throw new Error("invalid size");
        return dx = _0, dy = _1, resize();
    };
    
    density.cellSize = function(_) {
        if (!arguments.length) return 1 << k;
        if (!((_ = +_) >= 1)) throw new Error("invalid cell size");
        return k = Math.floor(Math.log(_) / Math.LN2), resize();
    };
    
    density.thresholds = function(_) {
        return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? () => _.slice() : () => _, density) : threshold;
    };
    
    density.bandwidth = function(_) {
        if (!arguments.length) return Math.sqrt(r * (1 << k));
        if (!((_ = +_) >= 0)) throw new Error("invalid bandwidth");
        r = _ * _ / (1 << k);
        return density;
    };
    
    density.mode = function(_) {
        if (!arguments.length) return mode;
        if (!["lines", "surfaces", "both"].includes(_)) throw new Error("invalid mode");
        mode = _;
        return density;
    };
    
    density.nullValuesMask = function(_) {
        if (!arguments.length) return nullValuesMask;
        nullValuesMask = !!_;
        return density;
    };
    
    density.bandwidthAdjust = function(_) {
        if (!arguments.length) return bandwidthAdjust;
        bandwidthAdjust = +_;
        if (bandwidthAdjust <= 0) bandwidthAdjust = 1.0;
        return density;
    };
    
    function resize() {
        o = r * 3;
        n = (dx + o * 2) >> k;
        m = (dy + o * 2) >> k;
        return density;
    }
    
    return density;
}

exports.contourDensity = density;
exports.contours = contours;
exports.toGeoJSON = toGeoJSON;

}));
