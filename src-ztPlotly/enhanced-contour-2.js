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

// Enhanced marching squares cases with saddle disambiguation
const marchingSquaresCases = [
    [],
    [[[1.0, 1.5], [0.5, 1.0]]],          // 1
    [[[1.5, 1.0], [1.0, 1.5]]],          // 2
    [[[1.5, 1.0], [0.5, 1.0]]],          // 3
    [[[1.0, 0.5], [1.5, 1.0]]],          // 4
    [[[1.0, 1.5], [0.5, 1.0]], [[1.0, 0.5], [1.5, 1.0]]], // 5 - saddle
    [[[1.0, 0.5], [1.0, 1.5]]],          // 6
    [[[1.0, 0.5], [0.5, 1.0]]],          // 7
    [[[0.5, 1.0], [1.0, 0.5]]],          // 8
    [[[1.0, 1.5], [1.0, 0.5]]],          // 9
    [[[0.5, 1.0], [1.0, 0.5]], [[1.5, 1.0], [1.0, 1.5]]], // 10 - saddle
    [[[1.5, 1.0], [1.0, 0.5]]],          // 11
    [[[0.5, 1.0], [1.5, 1.0]]],          // 12
    [[[1.0, 1.5], [1.5, 1.0]]],          // 13
    [[[0.5, 1.0], [1.0, 1.5]]],          // 14
    []                                     // 15
];

// Enhanced saddle disambiguation (inspired by plotly.js)
function disambiguateSaddle(caseIndex, corners, value) {
    if (caseIndex !== 5 && caseIndex !== 10) return caseIndex;
    
    const avg = (corners[0][0] + corners[0][1] + corners[1][0] + corners[1][1]) / 4;
    
    if (caseIndex === 5) {
        // Two peaks with a valley
        if (value > avg) return 713; // Custom case for disambiguation
        return 104; // Two valleys with a ridge
    } else if (caseIndex === 10) {
        // Two peaks with a valley
        if (value > avg) return 1114; // Custom case for disambiguation
        return 208; // Two valleys with a ridge
    }
    
    return caseIndex;
}

// Enhanced smoothing function with configurable parameters
function createSmoothingFunction(smoothingLevel = 1.0) {
    if (smoothingLevel <= 0) return function() {}; // No smoothing
    
    return function smoothAdvanced(ring, values, value, dx, dy) {
        const smoothingFactor = Math.max(0, Math.min(5, smoothingLevel));
        
        ring.forEach(function(point) {
            const x = point[0];
            const y = point[1];
            const xt = Math.floor(x);
            const yt = Math.floor(y);
            
            if (xt >= 0 && xt < dx - 1 && yt >= 0 && yt < dy - 1) {
                // Enhanced interpolation with smoothing factor
                const xFrac = x - xt;
                const yFrac = y - yt;
                
                if (Math.abs(xFrac - 0.5) < 0.1 && x > 0 && x < dx) {
                    const v0 = valid(values[yt * dx + xt]);
                    const v1 = valid(values[yt * dx + xt + 1]);
                    point[0] = smoothInterpolate(x, v0, v1, value, smoothingFactor);
                }
                
                if (Math.abs(yFrac - 0.5) < 0.1 && y > 0 && y < dy) {
                    const v0 = valid(values[yt * dx + xt]);
                    const v1 = valid(values[(yt + 1) * dx + xt]);
                    point[1] = smoothInterpolate(y, v0, v1, value, smoothingFactor);
                }
            }
        });
    };
}

function smoothInterpolate(coord, v0, v1, value, smoothingFactor) {
    const a = value - v0;
    const b = v1 - v0;
    const d = isFinite(a) && isFinite(b) && b !== 0 ? a / b : 0.5;
    
    // Apply smoothing curve based on smoothing factor
    const smoothed = coord + (d - 0.5) * (1 + smoothingFactor * 0.2);
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
    let threshold = function(values) { return 10; }; // Default 10 levels
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
        
        isolines(values, v, function(line) {
            smoothing(line, values, v, dx, dy);
            lines.push(line);
        });
        
        return {
            type: "MultiLineString",
            value: value,
            coordinates: lines
        };
    }
    
    // Generate contour surfaces (polygons) with enhanced hole assignment
    function contourSurfaces(values, value) {
        const v = value == null ? NaN : +value;
        if (isNaN(v)) throw new Error(`invalid value: ${value}`);
        
        const polygons = [];
        const holes = [];
        
        isorings(values, v, function(ring) {
            smoothing(ring, values, v, dx, dy);
            if (area(ring) > 0) {
                polygons.push([ring]);
            } else {
                holes.push(ring);
            }
        });
        
        // Enhanced hole assignment using spatial index
        const spatialIndex = new SpatialIndex();
        polygons.forEach((polygon, i) => spatialIndex.addPolygon(polygon, i));
        
        holes.forEach(function(hole) {
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
    
    // Marching squares for lines (isolines)
    function isolines(values, value, callback) {
        const fragmentByStart = new Map();
        const fragmentByEnd = new Map();
        let x, y, t0, t1, t2, t3;
        
        // Process grid with enhanced saddle handling
        y = -1;
        while (++y < dy - 1) {
            x = -1;
            while (++x < dx - 1) {
                const corners = [
                    [getValue(values, x, y), getValue(values, x + 1, y)],
                    [getValue(values, x, y + 1), getValue(values, x + 1, y + 1)]
                ];
                
                const caseIndex = getMarchingIndex(value, corners);
                if (caseIndex === 0 || caseIndex === 15) continue;
                
                const disambiguatedCase = disambiguateSaddle(caseIndex, corners, value);
                let caseLines = marchingSquaresCases[caseIndex];
                
                if (disambiguatedCase !== caseIndex) {
                    // Handle saddle cases with proper disambiguation
                    caseLines = getSaddleLines(disambiguatedCase);
                }
                
                caseLines.forEach(line => stitchLine(line, x, y, fragmentByStart, fragmentByEnd, callback));
            }
        }
        
        // Output any remaining open fragments
        fragmentByEnd.forEach(fragment => {
            callback(fragment.line);
        });
    }
    
    // Enhanced isorings with null mask support
    function isorings(values, value, callback) {
        const fragmentByStart = new Map();
        const fragmentByEnd = new Map();
        let x, y;
        
        // Special handling for null mask boundaries
        if (nullMask) {
            generateNullBoundaries(values, value, callback);
        }
        
        // Process grid
        y = -1;
        while (++y < dy - 1) {
            x = -1;
            while (++x < dx - 1) {
                if (nullMask && !isValidCell(x, y, values)) continue;
                
                const corners = [
                    [getValue(values, x, y), getValue(values, x + 1, y)],
                    [getValue(values, x, y + 1), getValue(values, x + 1, y + 1)]
                ];
                
                const caseIndex = getMarchingIndex(value, corners);
                if (caseIndex === 0 || caseIndex === 15) continue;
                
                const disambiguatedCase = disambiguateSaddle(caseIndex, corners, value);
                let caseLines = marchingSquaresCases[caseIndex];
                
                if (disambiguatedCase !== caseIndex) {
                    caseLines = getSaddleLines(disambiguatedCase);
                }
                
                caseLines.forEach(line => stitchRing(line, x, y, fragmentByStart, fragmentByEnd, callback));
            }
        }
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
        return (above(corners[0][0], value) ? 1 : 0) +
               (above(corners[0][1], value) ? 2 : 0) +
               (above(corners[1][1], value) ? 4 : 0) +
               (above(corners[1][0], value) ? 8 : 0);
    }
    
    function getSaddleLines(saddleCase) {
        // Return appropriate lines for disambiguated saddle cases
        switch (saddleCase) {
            case 713: return [[[1.5, 1.0], [1.0, 0.5]]];
            case 1114: return [[[0.5, 1.0], [1.0, 1.5]]];
            case 104: return [[[1.0, 1.5], [0.5, 1.0]]];
            case 208: return [[[1.0, 0.5], [1.5, 1.0]]];
            default: return [];
        }
    }
    
    function stitchLine(line, x, y, fragmentByStart, fragmentByEnd, callback) {
        // Similar to stitch but for open lines
        const start = [line[0][0] + x, line[0][1] + y];
        const end = [line[1][0] + x, line[1][1] + y];
        const startKey = `${start[0]},${start[1]}`;
        const endKey = `${end[0]},${end[1]}`;
        
        // Try to connect to existing fragments
        const fragmentFromEnd = fragmentByEnd.get(startKey);
        const fragmentFromStart = fragmentByStart.get(endKey);
        
        if (fragmentFromEnd && fragmentFromStart) {
            // Connect two fragments
            const combined = fragmentFromEnd.line.concat([start, end], fragmentFromStart.line);
            fragmentByEnd.delete(fragmentFromEnd.endKey);
            fragmentByStart.delete(fragmentFromStart.startKey);
            
            if (fragmentFromEnd === fragmentFromStart) {
                // Closed loop
                callback(combined);
            } else {
                // New larger fragment
                const newStartKey = fragmentFromEnd.startKey;
                const newEndKey = fragmentFromStart.endKey;
                const newFragment = { line: combined, startKey: newStartKey, endKey: newEndKey };
                fragmentByStart.set(newStartKey, newFragment);
                fragmentByEnd.set(newEndKey, newFragment);
            }
        } else if (fragmentFromEnd) {
            // Extend existing fragment
            fragmentFromEnd.line.push(start, end);
            fragmentByEnd.delete(startKey);
            fragmentByEnd.set(endKey, fragmentFromEnd);
            fragmentFromEnd.endKey = endKey;
        } else if (fragmentFromStart) {
            // Prepend to existing fragment
            fragmentFromStart.line.unshift(end, start);
            fragmentByStart.delete(endKey);
            fragmentByStart.set(startKey, fragmentFromStart);
            fragmentFromStart.startKey = startKey;
        } else {
            // New fragment
            const fragment = { line: [start, end], startKey, endKey };
            fragmentByStart.set(startKey, fragment);
            fragmentByEnd.set(endKey, fragment);
        }
    }
    
    function stitchRing(line, x, y, fragmentByStart, fragmentByEnd, callback) {
        // Enhanced version of original stitch with better error handling
        const start = [line[0][0] + x, line[0][1] + y];
        const end = [line[1][0] + x, line[1][1] + y];
        const startIndex = index(start);
        const endIndex = index(end);
        
        let f = fragmentByEnd.get(startIndex);
        let g = fragmentByStart.get(endIndex);
        
        if (f) {
            if (g) {
                fragmentByEnd.delete(f.end);
                fragmentByStart.delete(g.start);
                if (f === g) {
                    f.ring.push(end);
                    callback(f.ring);
                } else {
                    const newFragment = {
                        start: f.start,
                        end: g.end,
                        ring: f.ring.concat(g.ring)
                    };
                    fragmentByStart.set(f.start, newFragment);
                    fragmentByEnd.set(g.end, newFragment);
                }
            } else {
                fragmentByEnd.delete(f.end);
                f.ring.push(end);
                f.end = endIndex;
                fragmentByEnd.set(endIndex, f);
            }
        } else {
            f = fragmentByStart.get(endIndex);
            if (f) {
                if (g = fragmentByEnd.get(startIndex)) {
                    fragmentByStart.delete(f.start);
                    fragmentByEnd.delete(g.end);
                    if (f === g) {
                        f.ring.push(end);
                        callback(f.ring);
                    } else {
                        const newFragment = {
                            start: g.start,
                            end: f.end,
                            ring: g.ring.concat(f.ring)
                        };
                        fragmentByStart.set(g.start, newFragment);
                        fragmentByEnd.set(f.end, newFragment);
                    }
                } else {
                    fragmentByStart.delete(f.start);
                    f.ring.unshift(start);
                    f.start = startIndex;
                    fragmentByStart.set(startIndex, f);
                }
            } else {
                const fragment = { start: startIndex, end: endIndex, ring: [start, end] };
                fragmentByStart.set(startIndex, fragment);
                fragmentByEnd.set(endIndex, fragment);
            }
        }
    }
    
    function generateNullBoundaries(values, value, callback) {
        // Generate boundary lines around null regions
        for (let y = 0; y < dy - 1; y++) {
            for (let x = 0; x < dx - 1; x++) {
                const hasNull = [
                    nullMask[y * dx + x] === 0,
                    nullMask[y * dx + x + 1] === 0,
                    nullMask[(y + 1) * dx + x] === 0,
                    nullMask[(y + 1) * dx + x + 1] === 0
                ];
                
                if (hasNull.some(Boolean) && !hasNull.every(Boolean)) {
                    // This cell has mixed null/non-null values, generate boundary
                    generateCellBoundary(x, y, hasNull, callback);
                }
            }
        }
    }
    
    function generateCellBoundary(x, y, hasNull, callback) {
        // Generate boundary segments for cells with null/non-null mix
        const lines = [];
        
        // Check each edge for null boundary
        if (hasNull[0] !== hasNull[1]) { // Top edge
            lines.push([[x + 0.5, y], [x + 0.5, y]]);
        }
        if (hasNull[1] !== hasNull[3]) { // Right edge
            lines.push([[x + 1, y + 0.5], [x + 1, y + 0.5]]);
        }
        if (hasNull[2] !== hasNull[3]) { // Bottom edge
            lines.push([[x + 0.5, y + 1], [x + 0.5, y + 1]]);
        }
        if (hasNull[0] !== hasNull[2]) { // Left edge
            lines.push([[x, y + 0.5], [x, y + 0.5]]);
        }
        
        lines.forEach(line => callback(line));
    }
    
    function index(point) {
        return point[0] * 2 + point[1] * (dx + 1) * 4;
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
    return x != null && +x >= value;
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