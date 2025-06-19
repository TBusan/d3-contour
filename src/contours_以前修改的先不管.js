import {extent, nice, thresholdSturges, ticks} from "d3-array";
import {slice} from "./array.js";
import ascending from "./ascending.js";
import area from "./area.js";
import constant from "./constant.js";
import contains from "./contains.js";
import noop from "./noop.js";

// Define constants for saddle point handling, inspired by plotly.js approach
const CHOOSESADDLE = {
  713: [7, 13],
  1114: [11, 14],
  104: [1, 4],
  208: [2, 8]
};
const SADDLEREMAINDER = {
  1: 4,
  2: 8,
  4: 1,
  7: 13,
  8: 2,
  11: 14,
  13: 7,
  14: 11
};

// The marching squares algorithm produces square-based contour rings with 16 cases
var cases = [
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

export default function() {
  var dx = 1,
      dy = 1,
      threshold = thresholdSturges,
      smooth = smoothLinear;

  function contours(values) {
    var tz = threshold(values);

    // Convert number of thresholds into uniform thresholds.
    if (!Array.isArray(tz)) {
      const e = extent(values, finite);
      tz = ticks(...nice(e[0], e[1], tz), tz);
      while (tz[tz.length - 1] >= e[1]) tz.pop();
      while (tz[1] < e[0]) tz.shift();
    } else {
      tz = tz.slice().sort(ascending);
    }

    return tz.map(value => contour(values, value));
  }

  // Accumulate, smooth contour rings, assign holes to exterior rings.
  // Based on https://github.com/mbostock/shapefile/blob/v0.6.2/shp/polygon.js
  function contour(values, value) {
    const v = value == null ? NaN : +value;
    if (isNaN(v)) throw new Error(`invalid value: ${value}`);

    var polygons = [],
        holes = [];

    isorings(values, v, function(ring) {
      smooth(ring, values, v);
      // Add distance information for post-processing
      addGridInfo(ring);
      if (area(ring) > 0) polygons.push([ring]);
      else holes.push(ring);
    });

    // Post-process rings to improve quality
    polygons = polygons.map(polygon => {
      // Remove too close points
      polygon[0] = optimizePointDensity(polygon[0]);
      return polygon;
    });

    holes = holes.map(hole => optimizePointDensity(hole));

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

  // Add grid index information to ring points for better path optimization
  function addGridInfo(ring) {
    ring.forEach(function(point) {
      point[2] = point[0] | 0; // add grid x index
      point[3] = point[1] | 0; // add grid y index
    });
  }

  // Remove points that are too close together
  function optimizePointDensity(ring) {
    if (ring.length <= 2) return ring;
    
    // Calculate total path length in grid units
    let totaldist = 0;
    const alldists = [];
    
    for (let i = 1; i < ring.length; i++) {
      const thisdist = ptDist(ring[i], ring[i - 1]);
      totaldist += thisdist;
      alldists.push(thisdist);
    }
    
    // Add last segment if this is a closed path
    const closedpath = equalPts(ring[0], ring[ring.length - 1], 0.01, 0.01);
    if (closedpath && ring.length > 1) {
      const thisdist = ptDist(ring[0], ring[ring.length - 1]);
      totaldist += thisdist;
      alldists.push(thisdist);
    }
    
    // If there are not enough points, return the original
    if (ring.length < 4) return ring;
    
    const distThresholdFactor = 0.2 * (smooth === smoothLinear ? 1 : 0);
    const distThreshold = totaldist / alldists.length * distThresholdFactor;
    
    // Skip optimization if threshold is zero or too small
    if (distThreshold <= 0.01) return ring;
    
    const result = [];
    let i = 0;
    let current = ring[0];
    result.push(current);
    
    while (i < ring.length - 1) {
      let nextIndex = i + 1;
      let distAcc = alldists[i];
      
      // Accumulate points that are too close
      while (nextIndex < ring.length - 1 && distAcc + alldists[nextIndex] < distThreshold) {
        distAcc += alldists[nextIndex];
        nextIndex++;
      }
      
      if (nextIndex > i + 1) {
        // Average the points
        const avgPoint = [0, 0];
        for (let j = i + 1; j <= nextIndex; j++) {
          avgPoint[0] += ring[j][0];
          avgPoint[1] += ring[j][1];
        }
        avgPoint[0] /= (nextIndex - i);
        avgPoint[1] /= (nextIndex - i);
        result.push(avgPoint);
      } else {
        result.push(ring[nextIndex]);
      }
      
      i = nextIndex;
    }
    
    // Ensure closed paths remain closed
    if (closedpath && !equalPts(result[0], result[result.length - 1], 0.01, 0.01)) {
      result.push([result[0][0], result[0][1]]);
    }
    
    return result;
  }

  // Marching squares with isolines stitched into rings.
  // Based on https://github.com/topojson/topojson-client/blob/v3.0.0/src/stitch.js
  // Enhanced with better saddle point handling
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
      t1 = above(values[y * dx + dx], value);
      t2 = above(values[y * dx], value);
      cases[t1 << 1 | t2 << 2].forEach(stitch);
      while (++x < dx - 1) {
        t0 = t1, t1 = above(values[y * dx + dx + x + 1], value);
        t3 = t2, t2 = above(values[y * dx + x + 1], value);
        
        // Get marching index
        let mi = t0 | t1 << 1 | t2 << 2 | t3 << 3;
        
        // Handle saddle cases better (case 5 and 10)
        if (mi === 5 || mi === 10) {
          // Sample the center value to disambiguate the saddle
          const corners = [
            [values[y * dx + x], values[y * dx + x + 1]],
            [values[(y + 1) * dx + x], values[(y + 1) * dx + x + 1]]
          ];
          const avg = (corners[0][0] + corners[0][1] + corners[1][0] + corners[1][1]) / 4;
          
          // Determine which way the saddle should go based on the average
          if (value > avg) {
            mi = (mi === 5) ? 713 : 1114;
          } else {
            mi = (mi === 5) ? 104 : 208;
          }
          
          // Handle the special saddle cases
          if (mi > 15) {
            if (mi === 713 || mi === 1114) {
              // Apply the saddle case direction using CHOOSESADDLE
              const directions = CHOOSESADDLE[mi];
              cases[directions[0]].forEach(stitch);
              cases[directions[1]].forEach(stitch);
              
              // Update the marching index for potential next steps
              mi = SADDLEREMAINDER[directions[1]];
            } else if (mi === 104 || mi === 208) {
              // Apply the saddle case direction using CHOOSESADDLE
              const directions = CHOOSESADDLE[mi];
              cases[directions[0]].forEach(stitch);
              cases[directions[1]].forEach(stitch);
              
              // Update the marching index for potential next steps
              mi = SADDLEREMAINDER[directions[1]];
            }
            continue;
          }
        }
        
        cases[mi].forEach(stitch);
      }
      cases[t1 | t2 << 3].forEach(stitch);
    }

    // Special case for the last row (y = dy - 1, t0 = t1 = 0).
    x = -1;
    t2 = values[y * dx] >= value;
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

  // Compare if two points are equal within tolerance
  function equalPts(pt1, pt2, xtol, ytol) {
    return Math.abs(pt1[0] - pt2[0]) < xtol &&
           Math.abs(pt1[1] - pt2[1]) < ytol;
  }

  // Calculate distance between points in grid units
  function ptDist(pt1, pt2) {
    // If grid indices are available, use them
    if (pt1.length > 2 && pt2.length > 2) {
      const dx = pt1[2] - pt2[2];
      const dy = pt1[3] - pt2[3];
      return Math.sqrt(dx * dx + dy * dy);
    } else {
      // Fall back to pixel coordinates
      const dx = pt1[0] - pt2[0];
      const dy = pt1[1] - pt2[1];
      return Math.sqrt(dx * dx + dy * dy);
    }
  }

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

// When computing the extent, ignore infinite values (as well as invalid ones).
function finite(x) {
  return isFinite(x) ? x : NaN;
}

// Is the (possibly invalid) x greater than or equal to the (known valid) value?
// Treat any invalid value as below negative infinity.
function above(x, value) {
  return x == null ? false : +x >= value;
}

// During smoothing, treat any invalid value as negative infinity.
function valid(v) {
  return v == null || isNaN(v = +v) ? -Infinity : v;
}

function smooth1(x, v0, v1, value) {
  const a = value - v0;
  const b = v1 - v0;
  const d = isFinite(a) || isFinite(b) ? a / b : Math.sign(a) / Math.sign(b);
  return isNaN(d) ? x : x + d - 0.5;
}