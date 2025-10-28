/**
 * Fixed Enhanced Contour Library
 * A more faithful enhancement of d3-contour with careful improvements
 * 
 * Features:
 * - Support for both contour lines and contour surfaces
 * - Proper saddle point disambiguation
 * - More accurate isocontour tracing
 * - Optimized hole assignment algorithm
 * - Null value masking
 */

import { extent, ticks, thresholdSturges } from "d3-array";

// Utility functions
function ascending(a, b) {
    return a - b;
}

// Enhanced marching squares cases - matched exactly to d3-contour coordinates
const cases = [
    [],
    [[[1.0, 1.5], [0.5, 1.0]]],
    [[[1.5, 1.0], [1.0, 1.5]]],
    [[[1.5, 1.0], [0.5, 1.0]]],
    [[[1.0, 0.5], [1.5, 1.0]]],
    [[[1.0, 1.5], [0.5, 1.0]], [[1.0, 0.5], [1.5, 1.0]]],
    [[[1.0, 0.5], [1.0, 1.5]]],
    [[[1.0, 0.5], [0.5, 1.0]]],
    [[[0.5, 1.0], [1.0, 0.5]]],
    [[[1.0, 1.5], [1.0, 0.5]]],
    [[[0.5, 1.0], [1.0, 0.5]], [[1.5, 1.0], [1.0, 1.5]]],
    [[[1.5, 1.0], [1.0, 0.5]]],
    [[[0.5, 1.0], [1.5, 1.0]]],
    [[[1.0, 1.5], [1.5, 1.0]]],
    [[[0.5, 1.0], [1.0, 1.5]]],
    []
];

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

// Core contour generator
export default function() {
    let dx = 1;
    let dy = 1;
    let threshold = thresholdSturges;
    let smooth = smoothLinear;
    let mode = 'surfaces'; // 'surfaces', 'lines', or 'both'
    let nullMask = null;
    
    function contours(values) {
        const tz = computeThresholds(values);
        
        if (mode === "lines") {
            return tz.map(value => contourLines(values, value));
        } else if (mode === "surfaces" || mode === undefined) {
            return tz.map(value => contour(values, value));
        } else { // both
            return tz.map(value => ({
                lines: contourLines(values, value),
                surfaces: contour(values, value),
                value: value
            }));
        }
    }
    
    function computeThresholds(values) {
        let tz = threshold(values);
        
        // Convert number of thresholds into uniform thresholds
        if (!Array.isArray(tz)) {
            const e = extent(values, finite);
            if (isNaN(e[0]) || isNaN(e[1])) return [];
            tz = ticks(e[0], e[1], tz);
        } else {
            tz = tz.slice().sort(ascending);
        }
        
        return tz;
    }
    
    // Original d3-contour algorithm for surfaces
    function contour(values, value) {
        const v = value == null ? NaN : +value;
        if (isNaN(v)) throw new Error(`invalid value: ${value}`);

        var polygons = [],
            holes = [];

        isorings(values, v, function(ring) {
            smooth(ring, values, v);
            if (area(ring) > 0) polygons.push([ring]);
            else holes.push(ring);
        });

        // Use optimized spatial index for hole assignment
        if (polygons.length > 0 && holes.length > 0) {
            const spatialIndex = new SpatialIndex();
            polygons.forEach((polygon, i) => spatialIndex.addPolygon(polygon, i));
            
            holes.forEach(function(hole) {
                const containingIndex = spatialIndex.findContainingPolygon(hole);
                if (containingIndex !== -1) {
                    polygons[containingIndex].push(hole);
                }
            });
        }

        return {
            type: "MultiPolygon",
            value: value,
            coordinates: polygons
        };
    }
    
    // Enhanced contour lines generator
    function contourLines(values, value) {
        const v = value == null ? NaN : +value;
        if (isNaN(v)) throw new Error(`invalid value: ${value}`);
        
        const lines = [];
        
        // Use the same isorings algorithm but collect open lines
        isoSegments(values, v, function(line) {
            smooth(line, values, v);
            lines.push(line);
        });
        
        return {
            type: "MultiLineString",
            value: value,
            coordinates: lines
        };
    }
    
    // Collect line segments without forming rings
    function isoSegments(values, value, callback) {
        var segments = [];
        
        // Trace through the grid using marching squares
        marchingSquares(values, value, (x, y, caseIndex) => {
            const caseSegments = cases[caseIndex];
            caseSegments.forEach(segment => {
                const line = segment.map(point => [
                    point[0] + x,
                    point[1] + y
                ]);
                segments.push(line);
            });
        });
        
        // Connect segments into continuous lines
        while (segments.length > 0) {
            const line = [segments[0][0], segments[0][1]];
            segments.splice(0, 1);
            
            // Try to extend the line
            let extended = true;
            while (extended) {
                extended = false;
                
                const end = line[line.length - 1];
                for (let i = 0; i < segments.length; i++) {
                    const segment = segments[i];
                    const start = segment[0];
                    const finish = segment[1];
                    
                    // If we can connect to the start of this segment
                    if (distance(end, start) < 1e-6) {
                        line.push(finish);
                        segments.splice(i, 1);
                        extended = true;
                        break;
                    }
                    // If we can connect to the end of this segment
                    else if (distance(end, finish) < 1e-6) {
                        line.push(start);
                        segments.splice(i, 1);
                        extended = true;
                        break;
                    }
                }
            }
            
            // Check if we found a closed loop
            if (distance(line[0], line[line.length - 1]) < 1e-6) {
                // For closed loops, make sure first and last points are identical
                line[line.length - 1] = line[0].slice();
            }
            
            callback(line);
        }
    }
    
    // Original d3-contour isorings algorithm
    function isorings(values, value, callback) {
        var fragmentByStart = new Array,
            fragmentByEnd = new Array;
        
        // Trace through the grid using marching squares
        marchingSquares(values, value, (x, y, caseIndex) => {
            const caseSegments = cases[caseIndex];
            caseSegments.forEach(segment => {
                const start = [segment[0][0] + x, segment[0][1] + y];
                const end = [segment[1][0] + x, segment[1][1] + y];
                
                stitch(start, end);
            });
        });
        
        function stitch(start, end) {
            var startIndex = index(start),
                endIndex = index(end),
                f, g;
                
            if (f = fragmentByEnd[startIndex]) {
                if (g = fragmentByStart[endIndex]) {
                    delete fragmentByEnd[f.end];
                    delete fragmentByStart[g.start];
                    if (f === g) {
                        f.ring.push(end);
                        callback(f.ring);
                    } else {
                        fragmentByStart[f.start] = fragmentByEnd[g.end] = {
                            start: f.start,
                            end: g.end,
                            ring: f.ring.concat(g.ring)
                        };
                    }
                } else {
                    delete fragmentByEnd[f.end];
                    f.ring.push(end);
                    fragmentByEnd[f.end = endIndex] = f;
                }
            } else if (f = fragmentByStart[endIndex]) {
                if (g = fragmentByEnd[startIndex]) {
                    delete fragmentByStart[f.start];
                    delete fragmentByEnd[g.end];
                    if (f === g) {
                        f.ring.push(end);
                        callback(f.ring);
                    } else {
                        fragmentByStart[g.start] = fragmentByEnd[f.end] = {
                            start: g.start,
                            end: f.end,
                            ring: g.ring.concat(f.ring)
                        };
                    }
                } else {
                    delete fragmentByStart[f.start];
                    f.ring.unshift(start);
                    fragmentByStart[f.start = startIndex] = f;
                }
            } else {
                fragmentByStart[startIndex] = fragmentByEnd[endIndex] = {
                    start: startIndex,
                    end: endIndex,
                    ring: [start, end]
                };
            }
        }
    }
    
    // Helper function for indexing points
    function index(point) {
        return point[0] * 2 + point[1] * (dx + 1) * 4;
    }
    
    // Core marching squares algorithm with saddle point disambiguation
    function marchingSquares(values, value, callback) {
        var x, y, t0, t1, t2, t3;

        // Special case for the first row (y = -1, t2 = t3 = 0).
        x = y = -1;
        t1 = above(values[0], value);
        cases[t1 << 1].forEach(segment => callback(x, y, t1 << 1));
        while (++x < dx - 1) {
            t0 = t1, t1 = above(values[x + 1], value);
            cases[t0 | t1 << 1].forEach(segment => callback(x, y, t0 | t1 << 1));
        }
        cases[t1 << 0].forEach(segment => callback(x, y, t1 << 0));

        // General case for the intermediate rows.
        while (++y < dy - 1) {
            x = -1;
            t1 = above(values[y * dx + dx], value);
            t2 = above(values[y * dx], value);
            cases[t1 << 1 | t2 << 2].forEach(segment => callback(x, y, t1 << 1 | t2 << 2));
            while (++x < dx - 1) {
                t0 = t1, t1 = above(values[y * dx + dx + x + 1], value);
                t3 = t2, t2 = above(values[y * dx + x + 1], value);
                
                let caseIndex = t0 | t1 << 1 | t2 << 2 | t3 << 3;
                
                // Saddle point disambiguation
                if (caseIndex === 5 || caseIndex === 10) {
                    // Calculate average value at center of the cell
                    const avg = (
                        getValue(values, x, y) + 
                        getValue(values, x + 1, y) + 
                        getValue(values, x, y + 1) + 
                        getValue(values, x + 1, y + 1)
                    ) / 4;
                    
                    // Select case based on center value comparison
                    if (caseIndex === 5) {
                        caseIndex = value > avg ? 5 : 10;
                    } else { // caseIndex === 10
                        caseIndex = value > avg ? 10 : 5;
                    }
                }
                
                cases[caseIndex].forEach(segment => callback(x, y, caseIndex));
            }
            cases[t1 | t2 << 3].forEach(segment => callback(x, y, t1 | t2 << 3));
        }

        // Special case for the last row (y = dy - 1, t0 = t1 = 0).
        x = -1;
        t2 = above(values[y * dx], value);
        cases[t2 << 2].forEach(segment => callback(x, y, t2 << 2));
        while (++x < dx - 1) {
            t3 = t2, t2 = above(values[y * dx + x + 1], value);
            cases[t2 << 2 | t3 << 3].forEach(segment => callback(x, y, t2 << 2 | t3 << 3));
        }
        cases[t2 << 3].forEach(segment => callback(x, y, t2 << 3));
    }
    
    function getValue(values, x, y) {
        if (x < 0 || x >= dx || y < 0 || y >= dy) return NaN;
        const index = y * dx + x;
        if (index < 0 || index >= values.length) return NaN;
        const value = values[index];
        return nullMask && nullMask[index] === 0 ? NaN : value;
    }
    
    // Original d3-contour smoothing function
    function smoothLinear(ring, values, value) {
        ring.forEach(function(point) {
            var x = point[0],
                y = point[1],
                xt = x | 0,
                yt = y | 0,
                v1 = valid(values[yt * dx + xt]);
            if (x > 0 && x < dx && xt === x) {
                point[0] = smooth1(x, valid(values[yt * dx + xt - 1]), v1, value);
            }
            if (y > 0 && y < dy && yt === y) {
                point[1] = smooth1(y, valid(values[(yt - 1) * dx + xt]), v1, value);
            }
        });
    }
    
    // API methods
    contours.contour = contour;
    
    contours.size = function(_) {
        if (!arguments.length) return [dx, dy];
        var _0 = Math.floor(_[0]), _1 = Math.floor(_[1]);
        if (!(_0 >= 0 && _1 >= 0)) throw new Error("invalid size");
        return dx = _0, dy = _1, contours;
    };
    
    contours.thresholds = function(_) {
        return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? () => _.slice() : () => _, contours) : threshold;
    };
    
    contours.smooth = function(_) {
        return arguments.length ? (smooth = _ ? smoothLinear : noop, contours) : smooth === smoothLinear;
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

// Helper functions
function noop() {}

function finite(x) {
    return isFinite(x) ? x : NaN;
}

function above(x, value) {
    return x != null && isFinite(x) && x >= value;
}

function valid(v) {
    return v == null || isNaN(v = +v) ? -Infinity : v;
}

function smooth1(x, v0, v1, value) {
    const a = value - v0;
    const b = v1 - v0;
    const d = isFinite(a) || isFinite(b) ? a / b : Math.sign(a) / Math.sign(b);
    return isNaN(d) ? x : x + d - 0.5;
}

function distance(p1, p2) {
    return Math.sqrt(Math.pow(p1[0] - p2[0], 2) + Math.pow(p1[1] - p2[1], 2));
}

// Area calculation (positive for counter-clockwise, negative for clockwise)
function area(ring) {
    let i = 0, n = ring.length;
    let area = ring[n - 1][1] * ring[0][0] - ring[n - 1][0] * ring[0][1];
    while (++i < n) area += ring[i - 1][1] * ring[i][0] - ring[i - 1][0] * ring[i][1];
    return area;
}

// Point in polygon test (for hole assignment)
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

// Utility function to convert contour output to GeoJSON
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