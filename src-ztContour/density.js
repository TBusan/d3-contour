/**
 * Enhanced contour density generator
 * Improves the original d3-contour density with better handling of smoothing and nulls
 */

import { blur2, max, ticks } from "d3-array";
import contours from "./contours.js";
import { transformContours } from "./utils/geoutils.js";

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
export default function() {
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
        blur2({data: values, width: n, height: m}, r * pow2k * bandwidthAdjust);
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
            const maxVal = max(values) / pow4k;
            if (maxVal <= 0) return []; // No data
            
            // Generate thresholds
            tz = ticks(Number.MIN_VALUE, maxVal, tz);
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
            get: () => max(values) / pow4k
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