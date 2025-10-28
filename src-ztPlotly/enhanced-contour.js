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
        // IMPROVED: Apply more aggressive smoothing for grid artifacts
        const SMOOTH_THRESHOLD = 0.05; // Increased from 0.1
        
        ring.forEach(function(point) {
            const x = point[0];
            const y = point[1];
            const xt = Math.floor(x);
            const yt = Math.floor(y);
            
            if (xt >= 0 && xt < dx - 1 && yt >= 0 && yt < dy - 1) {
                // Enhanced interpolation
                const xFrac = x - xt;
                const yFrac = y - yt;
                
                // IMPROVED: Detect and smooth grid-aligned points
                if (Math.abs(xFrac) < SMOOTH_THRESHOLD || Math.abs(xFrac - 1) < SMOOTH_THRESHOLD) {
                    // Smoothly adjust x coordinate away from grid lines
                    point[0] += (xFrac < 0.5 ? 1 : -1) * SMOOTH_THRESHOLD * smoothingLevel;
                }
                
                if (Math.abs(yFrac) < SMOOTH_THRESHOLD || Math.abs(yFrac - 1) < SMOOTH_THRESHOLD) {
                    // Smoothly adjust y coordinate away from grid lines
                    point[1] += (yFrac < 0.5 ? 1 : -1) * SMOOTH_THRESHOLD * smoothingLevel;
                }
                
                if (Math.abs(xFrac - 0.5) < SMOOTH_THRESHOLD && x > 0 && x < dx) {
                    const v0 = valid(values[yt * dx + xt]);
                    const v1 = valid(values[yt * dx + xt + 1]);
                    point[0] = smoothInterpolate(x, v0, v1, value);
                }
                
                if (Math.abs(yFrac - 0.5) < SMOOTH_THRESHOLD && y > 0 && y < dy) {
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
                if (segment.length < 2) return; // Skip invalid segments
                
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
        let connectedLines = connectSegments(lines);
        
        // IMPROVED: Filter out very short lines and grid artifacts
        connectedLines = connectedLines.filter(line => {
            // Remove short lines (likely artifacts)
            if (line.length < 4) return false;
            
            // Calculate the path length
            let pathLength = 0;
            for (let i = 1; i < line.length; i++) {
                pathLength += distance(line[i-1], line[i]);
            }
            
            // Filter out lines that are too short
            if (pathLength < 2.0) return false;
            
            // Filter out grid-aligned artifacts
            let verticalSegments = 0;
            let horizontalSegments = 0;
            
            for (let i = 1; i < line.length; i++) {
                const dx = Math.abs(line[i][0] - line[i-1][0]);
                const dy = Math.abs(line[i][1] - line[i-1][1]);
                
                if (dx < 1e-4) verticalSegments++;
                if (dy < 1e-4) horizontalSegments++;
            }
            
            // If more than 50% of segments are strictly vertical or horizontal, it's likely a grid artifact
            const totalSegments = line.length - 1;
            if (verticalSegments > totalSegments * 0.5 || horizontalSegments > totalSegments * 0.5) {
                return false;
            }
            
            return true;
        });
        
        // IMPROVED: Try to close open loops that are nearly closed
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
        
        // Apply additional smoothing to reduce sharp angles
        connectedLines = smoothLines(connectedLines);
        
        // IMPROVED: Remove self-intersections
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
        
        // Assign holes to their containing polygons
        const spatialIndex = new SpatialIndex();
        polygons.forEach((polygon, i) => spatialIndex.addPolygon(polygon, i));
        
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
                
                // Skip cells with NaN values
                if (isNaN(corners[0][0]) || isNaN(corners[0][1]) || 
                    isNaN(corners[1][0]) || isNaN(corners[1][1])) {
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
                let bestMatchDistance = 1e-4; // Threshold for considering points as matching
                
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
                    bestMatchDistance = 1e-4;
                    
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
            if (line.length > 3 && distance(line[0], line[line.length - 1]) < 1e-4) {
                // Make sure the first and last points are exactly the same for a proper closed loop
                line[line.length - 1] = line[0].slice();
            }
            
            // Add the completed line to the result
            result.push(line);
        }
        
        return result;
    }
    
    // Helper function to calculate distance between two points
    function distance(p1, p2) {
        const dx = p1[0] - p2[0];
        const dy = p1[1] - p2[1];
        return Math.sqrt(dx * dx + dy * dy);
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
            if (segmentsIntersect(a, b, c, d)) {
                return true;
            }
        }
    }
    
    return false;
}

// Check if two line segments intersect
function segmentsIntersect(a, b, c, d) {
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

// Check if two points are equal (within a small epsilon)
function pointsEqual(p1, p2) {
    const epsilon = 1e-6;
    return Math.abs(p1[0] - p2[0]) < epsilon && Math.abs(p1[1] - p2[1]) < epsilon;
}

// IMPROVED: New function to smooth lines to remove sharp angles
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
        smoothed.push(line[line.length - 1]);
        
        return smoothed;
    });
}

// IMPROVED: New function to remove self-intersections in a line
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
            if (segmentsIntersect(
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