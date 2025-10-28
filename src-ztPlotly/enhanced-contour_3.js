/**
 * Enhanced Contour Library
 * Improved version of d3-contour with plotly.js optimizations
 * 
 * Features:
 * - Support for both contour lines and contour surfaces
 * - Enhanced smoothing with configurable parameters
 * - Proper saddle point disambiguation
 * - Optimized hole assignment algorithm
 * - Null value masking and boundary clipping
 * - GeoJSON export for both lines and polygons
 */

// Utility functions
function extent(values, accessor = x => x) {
    let min = Infinity, max = -Infinity;
    for (let i = 0; i < values.length; i++) {
        const v = accessor(values[i]);
        if (v != null && isFinite(v)) {
            if (v < min) min = v;
            if (v > max) max = v;
        }
    }
    return min === Infinity ? [NaN, NaN] : [min, max];
}

function ticks(start, stop, count) {
    const step = (stop - start) / Math.max(1, count - 1);
    const result = [];
    for (let i = 0; i < count; i++) {
        result.push(start + i * step);
    }
    return result;
}

function ascending(a, b) {
    return a - b;
}

// Enhanced marching squares cases
// Each case is an array of line segments
// Each segment connects two points on the edges of a unit square
// The points are specified as [x, y] coordinates where:
// - (0,0) is bottom-left corner
// - (1,0) is bottom-right corner
// - (0,1) is top-left corner
// - (1,1) is top-right corner
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

// Enhanced smoothing function with configurable parameters
function createSmoothingFunction(smoothingLevel = 1.0) {
    if (smoothingLevel <= 0) return function() {}; // No smoothing
    
    return function smoothAdvanced(ring, values, value, dx, dy) {
        ring.forEach(function(point) {
            const x = point[0];
            const y = point[1];
            const xt = Math.floor(x);
            const yt = Math.floor(y);
            
            if (xt >= 0 && xt < dx - 1 && yt >= 0 && yt < dy - 1) {
                // Enhanced interpolation
                const xFrac = x - xt;
                const yFrac = y - yt;
                
                if (Math.abs(xFrac - 0.5) < 0.1 && x > 0 && x < dx) {
                    const v0 = valid(values[yt * dx + xt]);
                    const v1 = valid(values[yt * dx + xt + 1]);
                    point[0] = smoothInterpolate(x, v0, v1, value);
                }
                
                if (Math.abs(yFrac - 0.5) < 0.1 && y > 0 && y < dy) {
                    const v0 = valid(values[yt * dx + xt]);
                    const v1 = valid(values[(yt + 1) * dx + xt]);
                    point[1] = smoothInterpolate(y, v0, v1, value);
                }
            }
        });
    };
}

function smoothInterpolate(coord, v0, v1, value) {
    const a = value - v0;
    const b = v1 - v0;
    const d = isFinite(a) && isFinite(b) && b !== 0 ? a / b : 0.5;
    
    // Apply simple linear interpolation
    const smoothed = Math.floor(coord) + d;
    return smoothed;
}

// Optimized spatial index for hole assignment
class SpatialIndex {
    constructor() {
        this.polygons = [];
        this.bounds = [];
    }
    
    addPolygon(polygon, index) {
        this.polygons.push({ polygon, index });
        this.bounds.push(this.computeBounds(polygon[0])); // Outer ring
    }
    
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
    
    findContainingPolygon(hole) {
        const holeBounds = this.computeBounds(hole);
        const candidates = [];
        
        // First pass: bounding box test
        for (let i = 0; i < this.bounds.length; i++) {
            const bounds = this.bounds[i];
            if (holeBounds.minX >= bounds.minX && holeBounds.maxX <= bounds.maxX &&
                holeBounds.minY >= bounds.minY && holeBounds.maxY <= bounds.maxY) {
                candidates.push(i);
            }
        }
        
        // Second pass: actual containment test
        for (const idx of candidates) {
            if (contains(this.polygons[idx].polygon[0], hole) !== -1) {
                return idx;
            }
        }
        
        return -1;
    }
}

// Enhanced contour generator
export default function enhancedContours() {
    let dx = 1;
    let dy = 1;
    let threshold = function() { return 10; }; // Default 10 levels
    let smoothing = createSmoothingFunction(1.0);
    let mode = 'surfaces'; // 'surfaces', 'lines', or 'both'
    let nullMask = null;
    
    function contours(values) {
        const tz = computeThresholds(values);
        
        if (mode === 'lines') {
            return tz.map(value => contourLines(values, value));
        } else if (mode === 'surfaces') {
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
        
        if (!Array.isArray(tz)) {
            const e = extent(values, finite);
            if (isNaN(e[0]) || isNaN(e[1])) return [];
            tz = ticks(e[0], e[1], tz);
        } else {
            tz = tz.slice().sort(ascending);
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
                const line = segment.map(point => [
                    point[0] + x,
                    point[1] + y
                ]);
                
                // Apply smoothing if enabled
                if (smoothing) {
                    smoothing(line, values, v, dx, dy);
                }
                
                lines.push(line);
            });
        });
        
        // Connect line segments into continuous lines
        const connectedLines = connectSegments(lines);
        
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
                if (smoothing) {
                    smoothing(line, values, v, dx, dy);
                }
                
                lines.push(line);
            });
        });
        
        // Connect line segments into closed rings
        const rings = connectSegmentsToRings(lines);
        
        // Separate rings into polygons and holes based on winding order
        const polygons = [];
        const holes = [];
        
        rings.forEach(ring => {
            if (area(ring) > 0) {
                // Counter-clockwise rings are polygon exteriors
                polygons.push([ring]);
            } else {
                // Clockwise rings are holes
                holes.push(ring);
            }
        });
        
        // Assign holes to their containing polygons
        const spatialIndex = new SpatialIndex();
        polygons.forEach((polygon, i) => spatialIndex.addPolygon(polygon, i));
        
        holes.forEach(hole => {
            const containingIndex = spatialIndex.findContainingPolygon(hole);
            if (containingIndex !== -1) {
                polygons[containingIndex].push(hole);
            }
        });
        
        return {
            type: "MultiPolygon",
            value: value,
            coordinates: polygons
        };
    }
    
    // Core marching squares algorithm
    function marchingSquares(values, value, callback) {
        // Process each cell in the grid
        for (let y = 0; y < dy - 1; y++) {
            for (let x = 0; x < dx - 1; x++) {
                // Skip cells with null values if mask is provided
                if (nullMask && !isValidCell(x, y, values)) continue;
                
                // Get the values at the four corners of the cell
                const corners = [
                    [getValue(values, x, y), getValue(values, x + 1, y)],
                    [getValue(values, x, y + 1), getValue(values, x + 1, y + 1)]
                ];
                
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
    
    // Get contour line segments for a given case
    function getContourSegments(caseIndex, corners, value) {
        // Get the basic segments for this case
        const segments = marchingSquaresCases[caseIndex];
        
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
        
        // Ensure result is within [0,1] range
        return isFinite(t) ? Math.max(0, Math.min(1, t)) : 0.5;
    }
    
    // Connect line segments into continuous lines
    function connectSegments(segments) {
        if (segments.length === 0) return [];
        
        const result = [];
        const used = new Set();
        
        // Start with any unused segment
        for (let i = 0; i < segments.length; i++) {
            if (used.has(i)) continue;
            
            used.add(i);
            let line = segments[i].slice();
            
            // Try to extend the line in both directions
            let extended;
            do {
                extended = false;
                
                // Try to extend at the end
                const end = line[line.length - 1];
                for (let j = 0; j < segments.length; j++) {
                    if (used.has(j)) continue;
                    
                    const segment = segments[j];
                    const start = segment[0];
                    const last = segment[segment.length - 1];
                    
                    // Check if this segment connects to our line
                    if (pointsEqual(end, start)) {
                        // Add all points except the first (which is a duplicate)
                        line.push(...segment.slice(1));
                        used.add(j);
                        extended = true;
                        break;
                    } else if (pointsEqual(end, last)) {
                        // Add all points except the last (which is a duplicate), in reverse order
                        line.push(...segment.slice(0, -1).reverse());
                        used.add(j);
                        extended = true;
                        break;
                    }
                }
                
                if (!extended) {
                    // Try to extend at the start
                    const start = line[0];
                    for (let j = 0; j < segments.length; j++) {
                        if (used.has(j)) continue;
                        
                        const segment = segments[j];
                        const segStart = segment[0];
                        const segEnd = segment[segment.length - 1];
                        
                        // Check if this segment connects to our line
                        if (pointsEqual(start, segEnd)) {
                            // Add all points except the last (which is a duplicate), at the beginning
                            line.unshift(...segment.slice(0, -1));
                            used.add(j);
                            extended = true;
                            break;
                        } else if (pointsEqual(start, segStart)) {
                            // Add all points except the first (which is a duplicate), in reverse order, at the beginning
                            line.unshift(...segment.slice(1).reverse());
                            used.add(j);
                            extended = true;
                            break;
                        }
                    }
                }
            } while (extended);
            
            // Check if the line forms a closed loop
            if (line.length > 2 && pointsEqual(line[0], line[line.length - 1])) {
                // Remove the duplicate end point for closed loops
                line.pop();
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
        
        // Try to close each line into a ring
        lines.forEach(line => {
            const start = line[0];
            const end = line[line.length - 1];
            
            // Check if the line is already a ring
            if (pointsEqual(start, end)) {
                // It's already a ring, just remove the duplicate end point
                rings.push(line.slice(0, -1));
                return;
            }
            
            // Try to find a line that can close this one
            for (let i = 0; i < lines.length; i++) {
                const otherLine = lines[i];
                if (otherLine === line) continue;
                
                const otherStart = otherLine[0];
                const otherEnd = otherLine[otherLine.length - 1];
                
                // Check if connecting these lines would form a ring
                if (pointsEqual(end, otherStart) && pointsEqual(otherEnd, start)) {
                    // Combine the lines to form a ring, removing duplicate points
                    const ring = [...line, ...otherLine.slice(1, -1)];
                    rings.push(ring);
                    return;
                }
            }
            
            // If we can't form a ring, just add the line as is
            // This shouldn't happen with proper contour data, but we handle it just in case
            rings.push(line);
        });
        
        return rings;
    }
    
    // Check if two points are equal (within a small epsilon)
    function pointsEqual(p1, p2) {
        const epsilon = 1e-6;
        return Math.abs(p1[0] - p2[0]) < epsilon && Math.abs(p1[1] - p2[1]) < epsilon;
    }
    
    function getValue(values, x, y) {
        if (x < 0 || x >= dx || y < 0 || y >= dy) return NaN;
        const index = y * dx + x;
        if (index < 0 || index >= values.length) return NaN;
        const value = values[index];
        return nullMask && nullMask[index] === 0 ? NaN : value;
    }
    
    function isValidCell(x, y, values) {
        const v1 = getValue(values, x, y);
        const v2 = getValue(values, x + 1, y);
        const v3 = getValue(values, x, y + 1);
        const v4 = getValue(values, x + 1, y + 1);
        return !isNaN(v1) && !isNaN(v2) && !isNaN(v3) && !isNaN(v4);
    }
    
    function getMarchingIndex(value, corners) {
        // Ensure proper case calculation by using explicit comparison
        const c00 = corners[0][0];
        const c10 = corners[0][1];
        const c11 = corners[1][1];
        const c01 = corners[1][0];
        
        // Calculate index based on which corners are above the threshold
        return (above(c00, value) ? 1 : 0) +
               (above(c10, value) ? 2 : 0) +
               (above(c11, value) ? 4 : 0) +
               (above(c01, value) ? 8 : 0);
    }
    
    // API
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
        if (!arguments.length) return smoothing !== noop;
        smoothing = typeof _ === 'number' ? createSmoothingFunction(_) : (_ ? createSmoothingFunction(1.0) : noop);
        return contours;
    };
    
    contours.mode = function(_) {
        if (!arguments.length) return mode;
        if (!['lines', 'surfaces', 'both'].includes(_)) throw new Error("invalid mode");
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

// Helper functions
function area(ring) {
    let i = 0, n = ring.length;
    let area = ring[n - 1][1] * ring[0][0] - ring[n - 1][0] * ring[0][1];
    while (++i < n) area += ring[i - 1][1] * ring[i][0] - ring[i - 1][0] * ring[i][1];
    return area;
}

function contains(ring, hole) {
    let i = -1, n = hole.length, c;
    while (++i < n) if (c = ringContains(ring, hole[i])) return c;
    return 0;
}

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

function segmentContains(a, b, c) {
    const i = +(a[0] === b[0]);
    return collinear(a, b, c) && within(a[i], c[i], b[i]);
}

function collinear(a, b, c) {
    return (b[0] - a[0]) * (c[1] - a[1]) === (c[0] - a[0]) * (b[1] - a[1]);
}

function within(p, q, r) {
    return p <= q && q <= r || r <= q && q <= p;
}

function finite(x) {
    return isFinite(x) ? x : NaN;
}

function above(x, value) {
    // More robust comparison to handle edge cases
    return x != null && isFinite(x) && x >= value;
}

function valid(v) {
    return v == null || isNaN(v = +v) ? -Infinity : v;
}

function noop() {}

// Export helper for GeoJSON conversion
export function toGeoJSON(contourResult, properties = {}) {
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