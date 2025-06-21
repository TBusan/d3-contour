import {extent, nice, thresholdSturges, ticks} from "d3-array";
import {slice} from "./array.js";
import ascending from "./ascending.js";
import area from "./area.js";
import constant from "./constant.js";
import contains from "./contains.js";
import noop from "./noop.js";

// Marching squares case lookup table
// These define how to connect points to form contour lines for each possible corner configuration
var cases = [
  [], // Case 0: All corners are below threshold
  [[[1.0, 1.5], [0.5, 1.0]]], // Case 1
  [[[1.5, 1.0], [1.0, 1.5]]], // Case 2
  [[[1.5, 1.0], [0.5, 1.0]]], // Case 3
  [[[1.0, 0.5], [1.5, 1.0]]], // Case 4
  [[[1.0, 1.5], [0.5, 1.0]], [[1.0, 0.5], [1.5, 1.0]]], // Case 5: Saddle point
  [[[1.0, 0.5], [1.0, 1.5]]], // Case 6
  [[[1.0, 0.5], [0.5, 1.0]]], // Case 7
  [[[0.5, 1.0], [1.0, 0.5]]], // Case 8
  [[[1.0, 1.5], [1.0, 0.5]]], // Case 9
  [[[0.5, 1.0], [1.0, 0.5]], [[1.5, 1.0], [1.0, 1.5]]], // Case 10: Saddle point
  [[[1.5, 1.0], [1.0, 0.5]]], // Case 11
  [[[0.5, 1.0], [1.5, 1.0]]], // Case 12
  [[[1.0, 1.5], [1.5, 1.0]]], // Case 13
  [[[0.5, 1.0], [1.0, 1.5]]], // Case 14
  [] // Case 15: All corners are above threshold
];

export default function() {
  var dx = 1, // Width of each grid cell
      dy = 1, // Height of each grid cell
      threshold = thresholdSturges,
      smooth = smoothLinear;

  function contours(values) {
    var tz = threshold(values);

    // Convert number of thresholds into uniform thresholds.
    if (!Array.isArray(tz)) {
      const validValues = values.filter(v => v != null && !isNaN(v));
      const e = extent(validValues, finite);
      if (e[0] === undefined || e[1] === undefined) {
        return []; // No valid data points
      }
      tz = ticks(...nice(e[0], e[1], tz), tz);
      while (tz[tz.length - 1] >= e[1]) tz.pop();
      while (tz[1] < e[0]) tz.shift();
    } else {
      tz = tz.slice().sort(ascending);
    }

    return tz.map(value => contour(values, value));
  }

  // Accumulate, smooth contour rings, assign holes to exterior rings.
  function contour(values, value) {
    const v = value == null ? NaN : +value;
    if (isNaN(v)) throw new Error(`invalid value: ${value}`);

    var polygons = [],
        holes = [];

    isorings(values, v, function(ring) {
      if (ring.length > 2) { // Ignore degenerate rings
        smooth(ring, values, v);
        if (area(ring) > 0) polygons.push([ring]);
        else holes.push(ring);
      }
    });

    holes.forEach(function(hole) {
      for (var i = 0, n = polygons.length, polygon; i < n; ++i) {
        if (contains((polygon = polygons[i])[0], hole) !== -1) {
          polygon.push(hole);
          return;
        }
      }
    });

    return {
      type: "MultiPolygon",
      value: value,
      coordinates: polygons
    };
  }

  // Marching squares with isolines stitched into rings.
  function isorings(values, value, callback) {
    var fragmentByStart = new Array,
        fragmentByEnd = new Array,
        x, y, t0, t1, t2, t3;

    // Special case for the first row (y = -1, t2 = t3 = 0).
    x = y = -1;
    t1 = above(values[0], value);
    cases[t1 << 1].forEach(stitch);
    while (++x < dx - 1) {
      t0 = t1, t1 = above(values[x + 1], value);
      cases[t0 | t1 << 1].forEach(stitch);
    }
    cases[t1 << 0].forEach(stitch);

    // General case for the intermediate rows.
    while (++y < dy - 1) {
      x = -1;
      t1 = above(values[(y + 1) * dx], value);
      t2 = above(values[y * dx], value);
      cases[t1 << 1 | t2 << 2].forEach(stitch);
      while (++x < dx - 1) {
        t0 = t1, t1 = above(values[(y + 1) * dx + x + 1], value);
        t3 = t2, t2 = above(values[y * dx + x + 1], value);
        
        // Special handling for ambiguous cases (saddle points)
        const config = t0 | t1 << 1 | t2 << 2 | t3 << 3;
        if (config === 5 || config === 10) {
          // For saddle points, use the average of cell values to determine configuration
          const avg = averageValidValues([
            values[y * dx + x],
            values[y * dx + x + 1], 
            values[(y + 1) * dx + x], 
            values[(y + 1) * dx + x + 1]
          ]);
          
          // Use the average value compared to threshold to resolve ambiguity
          if ((config === 5 && avg < value) || (config === 10 && avg < value)) {
            // Original configuration
            cases[config].forEach(stitch);
          } else {
            // Alternate configuration
            cases[15 - config].forEach(stitch);
          }
        } else {
          // Normal cases
          cases[config].forEach(stitch);
        }
      }
      cases[t1 | t2 << 3].forEach(stitch);
    }

    // Special case for the last row (y = dy - 1, t0 = t1 = 0).
    x = -1;
    t2 = above(values[y * dx], value);
    cases[t2 << 2].forEach(stitch);
    while (++x < dx - 1) {
      t3 = t2, t2 = above(values[y * dx + x + 1], value);
      cases[t2 << 2 | t3 << 3].forEach(stitch);
    }
    cases[t2 << 3].forEach(stitch);

    function stitch(line) {
      var start = [line[0][0] + x, line[0][1] + y],
          end = [line[1][0] + x, line[1][1] + y],
          startIndex = index(start),
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
            fragmentByStart[f.start] = fragmentByEnd[g.end] = {start: f.start, end: g.end, ring: f.ring.concat(g.ring)};
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
            fragmentByStart[g.start] = fragmentByEnd[f.end] = {start: g.start, end: f.end, ring: g.ring.concat(f.ring)};
          }
        } else {
          delete fragmentByStart[f.start];
          f.ring.unshift(start);
          fragmentByStart[f.start = startIndex] = f;
        }
      } else {
        fragmentByStart[startIndex] = fragmentByEnd[endIndex] = {start: startIndex, end: endIndex, ring: [start, end]};
      }
    }
  }

  function index(point) {
    return point[0] * 2 + point[1] * (dx + 1) * 4;
  }

  function averageValidValues(values) {
    const validValues = values.filter(v => v != null && !isNaN(+v));
    if (validValues.length === 0) return -Infinity; // Default below threshold if no valid values
    return validValues.reduce((sum, v) => sum + (+v), 0) / validValues.length;
  }

  // Enhanced smoothing that handles null values
  function smoothLinear(ring, values, value) {
    ring.forEach(function(point) {
      var x = point[0],
          y = point[1],
          xt = x | 0,
          yt = y | 0;
          
      // Ensure coordinates are within bounds
      if (xt < 0 || yt < 0 || xt >= dx || yt >= dy) return;
      
      const idx = yt * dx + xt;
      const v1 = valid(values[idx]);
      
      if (x > 0 && x < dx && xt === x) {
        const leftIdx = yt * dx + (xt - 1);
        const leftVal = valid(values[leftIdx]);
        
        if (isFinite(leftVal) && isFinite(v1)) {
          point[0] = smooth1(x, leftVal, v1, value);
        }
      }
      
      if (y > 0 && y < dy && yt === y) {
        const topIdx = (yt - 1) * dx + xt;
        const topVal = valid(values[topIdx]);
        
        if (isFinite(topVal) && isFinite(v1)) {
          point[1] = smooth1(y, topVal, v1, value);
        }
      }
    });
  }

  contours.contour = contour;

  contours.size = function(_) {
    if (!arguments.length) return [dx, dy];
    var _0 = Math.floor(_[0]), _1 = Math.floor(_[1]);
    if (!(_0 >= 0 && _1 >= 0)) throw new Error("invalid size");
    return dx = _0, dy = _1, contours;
  };

  contours.thresholds = function(_) {
    return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? constant(slice.call(_)) : constant(_), contours) : threshold;
  };

  contours.smooth = function(_) {
    return arguments.length ? (smooth = _ ? smoothLinear : noop, contours) : smooth === smoothLinear;
  };

  return contours;
}

// When computing the extent, ignore infinite, null, and NaN values
function finite(x) {
  return isFinite(x) ? x : NaN;
}

// Is the value above or equal to threshold?
// Handle null, undefined, and NaN as below threshold
function above(x, value) {
  return x != null && !isNaN(x = +x) ? x >= value : false;
}

// During smoothing, treat any invalid value as negative infinity
function valid(v) {
  return v == null || isNaN(v = +v) ? -Infinity : v;
}

// Linear interpolation for smoothing
function smooth1(x, v0, v1, value) {
  const a = value - v0;
  const b = v1 - v0;
  const d = isFinite(a) || isFinite(b) ? a / b : Math.sign(a) / Math.sign(b);
  return isNaN(d) ? x : x + d - 0.5;
} 