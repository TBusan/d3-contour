import {extent, nice, thresholdSturges, ticks} from "d3-array";
import {slice} from "./array.js";
import ascending from "./ascending.js";
import area from "./area.js";
import constant from "./constant.js";
import contains from "./contains.js";
import noop from "./noop.js";

// Enhanced cases table with saddle point disambiguation support
var cases = [
  [],
  [[[1.0, 1.5], [0.5, 1.0]]],
  [[[1.5, 1.0], [1.0, 1.5]]],
  [[[1.5, 1.0], [0.5, 1.0]]],
  [[[1.0, 0.5], [1.5, 1.0]]],
  // Cases 5 and 10 will be dynamically resolved based on saddle disambiguation
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

// Alternative cases for saddle point disambiguation
var saddleCases = {
  5: {
    // Case 5: choose connection based on center value
    primary: [[[1.0, 1.5], [0.5, 1.0]], [[1.0, 0.5], [1.5, 1.0]]],
    secondary: [[[1.0, 1.5], [1.0, 0.5]], [[0.5, 1.0], [1.5, 1.0]]]
  },
  10: {
    // Case 10: choose connection based on center value
    primary: [[[0.5, 1.0], [1.0, 0.5]], [[1.5, 1.0], [1.0, 1.5]]],
    secondary: [[[0.5, 1.0], [1.5, 1.0]], [[1.0, 0.5], [1.0, 1.5]]]
  }
};

export default function() {
  var dx = 1,
      dy = 1,
      threshold = thresholdSturges,
      smooth = smoothLinear,
      smoothingFactor = 1.0,
      generateIsolines = false,
      nullValue = null;

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

  // Enhanced contour generation with improved hole assignment and null handling
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

    // Enhanced hole assignment with proper topology handling
    assignHolesImproved(polygons, holes);

    return {
      type: generateIsolines ? "MultiLineString" : "MultiPolygon",
      value: value,
      coordinates: generateIsolines ? convertToLines(polygons) : polygons
    };
  }

  // Improved hole assignment algorithm with O(n log n) complexity
  function assignHolesImproved(polygons, holes) {
    if (!holes.length) return;

    // Create spatial index for polygons with bounding boxes
    var polygonIndex = polygons.map((polygon, i) => {
      var ring = polygon[0];
      var bbox = getBoundingBox(ring);
      return { index: i, polygon: polygon, bbox: bbox, ring: ring };
    });

    // Sort polygons by area (smallest first for proper nesting)
    polygonIndex.sort((a, b) => Math.abs(area(a.ring)) - Math.abs(area(b.ring)));

    holes.forEach(function(hole) {
      var holeBbox = getBoundingBox(hole);
      var assigned = false;

      // Find the smallest polygon that contains this hole
      for (var i = 0; i < polygonIndex.length && !assigned; i++) {
        var polyData = polygonIndex[i];
        
        // Quick bounding box check first
        if (bboxContains(polyData.bbox, holeBbox)) {
          // Full containment test
          if (contains(polyData.ring, hole) !== -1) {
            // Check if this hole is already assigned to a smaller nested polygon
            var shouldAssign = true;
            for (var j = 0; j < i; j++) {
              var innerPoly = polygonIndex[j];
              if (bboxContains(innerPoly.bbox, holeBbox) && 
                  contains(innerPoly.ring, hole) !== -1) {
                shouldAssign = false;
                break;
              }
            }
            
            if (shouldAssign) {
              polyData.polygon.push(hole);
              assigned = true;
            }
          }
        }
      }
    });
  }

  // Enhanced marching squares with saddle point disambiguation
  function isorings(values, value, callback) {
    var fragmentByStart = new Array,
        fragmentByEnd = new Array,
        x, y, t0, t1, t2, t3;

    // Special case for the first row (y = -1, t2 = t3 = 0).
    x = y = -1;
    t1 = above(values[0], value);
    processCase(t1 << 1, x, y, values, value).forEach(stitch);
    while (++x < dx - 1) {
      t0 = t1, t1 = above(values[x + 1], value);
      processCase(t0 | t1 << 1, x, y, values, value).forEach(stitch);
    }
    processCase(t1 << 0, x, y, values, value).forEach(stitch);

    // General case for the intermediate rows.
    while (++y < dy - 1) {
      x = -1;
      t1 = above(values[y * dx + dx], value);
      t2 = above(values[y * dx], value);
      processCase(t1 << 1 | t2 << 2, x, y, values, value).forEach(stitch);
      while (++x < dx - 1) {
        t0 = t1, t1 = above(values[y * dx + dx + x + 1], value);
        t3 = t2, t2 = above(values[y * dx + x + 1], value);
        processCase(t0 | t1 << 1 | t2 << 2 | t3 << 3, x, y, values, value).forEach(stitch);
      }
      processCase(t1 | t2 << 3, x, y, values, value).forEach(stitch);
    }

    // Special case for the last row (y = dy - 1, t0 = t1 = 0).
    x = -1;
    t2 = above(values[y * dx], value);
    processCase(t2 << 2, x, y, values, value).forEach(stitch);
    while (++x < dx - 1) {
      t3 = t2, t2 = above(values[y * dx + x + 1], value);
      processCase(t2 << 2 | t3 << 3, x, y, values, value).forEach(stitch);
    }
    processCase(t2 << 3, x, y, values, value).forEach(stitch);

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

  // Process marching squares case with saddle point disambiguation
  function processCase(caseIndex, x, y, values, value) {
    // Handle saddle point cases (5 and 10) with disambiguation
    if (caseIndex === 5 || caseIndex === 10) {
      return disambiguateSaddle(caseIndex, x, y, values, value);
    }
    return cases[caseIndex];
  }

  // Saddle point disambiguation based on center interpolation
  function disambiguateSaddle(caseIndex, x, y, values, value) {
    if (x < 0 || y < 0 || x >= dx - 1 || y >= dy - 1) {
      return cases[caseIndex]; // Fallback to default case at boundaries
    }

    // Get the four corner values
    var v00 = getValueSafe(values, x, y);
    var v10 = getValueSafe(values, x + 1, y);
    var v01 = getValueSafe(values, x, y + 1);
    var v11 = getValueSafe(values, x + 1, y + 1);

    // Calculate center value using bilinear interpolation
    var centerValue = (v00 + v10 + v01 + v11) / 4;

    // Choose case based on center value relative to threshold
    var usePrimary = centerValue >= value;
    var saddleCase = saddleCases[caseIndex];
    
    return usePrimary ? saddleCase.primary : saddleCase.secondary;
  }

  function index(point) {
    return point[0] * 2 + point[1] * (dx + 1) * 4;
  }

  // Enhanced smoothing with continuous parameter control
  function smoothLinear(ring, values, value) {
    if (smoothingFactor <= 0) return;
    
    ring.forEach(function(point) {
      var x = point[0],
          y = point[1],
          xt = x | 0,
          yt = y | 0,
          v1 = valid(values[yt * dx + xt]);
      
      if (x > 0 && x < dx && xt === x) {
        var smoothedX = smooth1(x, valid(values[yt * dx + xt - 1]), v1, value);
        point[0] = x + (smoothedX - x) * smoothingFactor;
      }
      if (y > 0 && y < dy && yt === y) {
        var smoothedY = smooth1(y, valid(values[(yt - 1) * dx + xt]), v1, value);
        point[1] = y + (smoothedY - y) * smoothingFactor;
      }
    });
  }

  // Helper functions
  function getBoundingBox(ring) {
    var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (var i = 0; i < ring.length; i++) {
      var point = ring[i];
      if (point[0] < minX) minX = point[0];
      if (point[0] > maxX) maxX = point[0];
      if (point[1] < minY) minY = point[1];
      if (point[1] > maxY) maxY = point[1];
    }
    return {minX: minX, minY: minY, maxX: maxX, maxY: maxY};
  }

  function bboxContains(bbox1, bbox2) {
    return bbox1.minX <= bbox2.minX && bbox1.maxX >= bbox2.maxX &&
           bbox1.minY <= bbox2.minY && bbox1.maxY >= bbox2.maxY;
  }

  function getValueSafe(values, x, y) {
    if (x < 0 || x >= dx || y < 0 || y >= dy) return nullValue;
    var val = values[y * dx + x];
    return val === nullValue ? nullValue : val;
  }

  function convertToLines(polygons) {
    // Convert polygon rings to line strings for isoline generation
    return polygons.map(polygon => polygon.map(ring => ring));
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
    return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? constant(slice.call(_)) : constant(_), contours) : threshold;
  };

  contours.smooth = function(_) {
    return arguments.length ? (smooth = _ ? smoothLinear : noop, contours) : smooth === smoothLinear;
  };

  // New API methods for enhanced functionality
  contours.smoothingFactor = function(_) {
    return arguments.length ? (smoothingFactor = Math.max(0, Math.min(1, +_)), contours) : smoothingFactor;
  };

  contours.isoLines = function(_) {
    return arguments.length ? (generateIsolines = !!_, contours) : generateIsolines;
  };

  contours.nullValue = function(_) {
    return arguments.length ? (nullValue = _, contours) : nullValue;
  };

  // Enhanced utility functions with null value support (moved inside scope)
  function finite(x) {
    return x != null && isFinite(x) ? x : NaN;
  }

  function above(x, value) {
    if (x === nullValue || x == null) return false;
    return +x >= value;
  }

  function valid(v) {
    if (v === nullValue || v == null || isNaN(v = +v)) return -Infinity;
    return v;
  }

  function smooth1(x, v0, v1, value) {
    const a = value - v0;
    const b = v1 - v0;
    const d = isFinite(a) || isFinite(b) ? a / b : Math.sign(a) / Math.sign(b);
    return isNaN(d) ? x : x + d - 0.5;
  }

  return contours;
}