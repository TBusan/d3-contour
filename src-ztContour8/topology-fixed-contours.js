import {extent, nice, thresholdSturges, ticks} from "d3-array";
import {slice} from "./array.js";
import ascending from "./ascending.js";
import area from "./area.js";
import constant from "./constant.js";
import contains from "./contains.js";
import noop from "./noop.js";

// Enhanced Marching Squares cases with saddle point disambiguation
var cases = [
  [],                                    // 0: 0000
  [[[1.0, 1.5], [0.5, 1.0]]],          // 1: 0001
  [[[1.5, 1.0], [1.0, 1.5]]],          // 2: 0010
  [[[1.5, 1.0], [0.5, 1.0]]],          // 3: 0011
  [[[1.0, 0.5], [1.5, 1.0]]],          // 4: 0100
  [[[1.0, 1.5], [0.5, 1.0]], [[1.0, 0.5], [1.5, 1.0]]], // 5: 0101 - saddle case
  [[[1.0, 0.5], [1.0, 1.5]]],          // 6: 0110
  [[[1.0, 0.5], [0.5, 1.0]]],          // 7: 0111
  [[[0.5, 1.0], [1.0, 0.5]]],          // 8: 1000
  [[[1.0, 1.5], [1.0, 0.5]]],          // 9: 1001
  [[[0.5, 1.0], [1.0, 0.5]], [[1.5, 1.0], [1.0, 1.5]]], // 10: 1010 - saddle case
  [[[1.5, 1.0], [1.0, 0.5]]],          // 11: 1011
  [[[0.5, 1.0], [1.5, 1.0]]],          // 12: 1100
  [[[1.0, 1.5], [1.5, 1.0]]],          // 13: 1101
  [[[0.5, 1.0], [1.0, 1.5]]],          // 14: 1110
  []                                     // 15: 1111
];

// Alternative cases for saddle disambiguation
var saddleCases = {
  5: [
    [[[1.0, 1.5], [0.5, 1.0]], [[1.0, 0.5], [1.5, 1.0]]], // separate
    [[[1.0, 1.5], [1.0, 0.5], [1.5, 1.0], [0.5, 1.0]]]   // connected
  ],
  10: [
    [[[0.5, 1.0], [1.0, 0.5]], [[1.5, 1.0], [1.0, 1.5]]], // separate
    [[[0.5, 1.0], [1.5, 1.0], [1.0, 1.5], [1.0, 0.5]]]   // connected
  ]
};

export default function() {
  var dx = 1,
      dy = 1,
      threshold = thresholdSturges,
      smooth = smoothLinear,
      smoothFactor = 0.5,
      enableSaddleDisambiguation = true,
      exportGeoJSON = true,
      handleNullValues = true,
      preventOverlap = true,
      extendBoundaries = true;

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

    // Generate contours with topology fixing
    if (preventOverlap) {
      return generateNonOverlappingContours(values, tz);
    } else {
      return tz.map(value => contour(values, value));
    }
  }

  // Generate non-overlapping contours using layered approach
  function generateNonOverlappingContours(values, thresholds) {
    const results = [];
    
    // Extend boundaries if needed
    const extendedValues = extendBoundaries ? extendGridBoundaries(values) : values;
    const extendedDx = extendBoundaries ? dx + 2 : dx;
    const extendedDy = extendBoundaries ? dy + 2 : dy;
    
    // Create level assignment for each cell
    const levelMap = new Array(extendedDx * extendedDy);
    
    for (let i = 0; i < extendedDx * extendedDy; i++) {
      const value = extendedValues[i];
      if (value == null || !isFinite(value)) {
        levelMap[i] = -1; // Mark as invalid
        continue;
      }
      
      // Find the highest threshold this value exceeds
      let level = -1;
      for (let j = 0; j < thresholds.length; j++) {
        if (value >= thresholds[j]) {
          level = j;
        } else {
          break;
        }
      }
      levelMap[i] = level;
    }
    
    // Generate contours for each level boundary
    for (let level = 0; level < thresholds.length; level++) {
      const contour = generateLevelContour(extendedValues, levelMap, thresholds[level], level, extendedDx, extendedDy);
      if (contour.coordinates.length > 0) {
        results.push(contour);
      }
    }
    
    return results;
  }

  // Generate contour for a specific level
  function generateLevelContour(values, levelMap, threshold, level, gridDx, gridDy) {
    var polygons = [],
        holes = [];

    // Temporarily update grid dimensions for extended grid
    const originalDx = dx, originalDy = dy;
    dx = gridDx; dy = gridDy;

    isoringsForLevel(values, levelMap, threshold, level, function(ring) {
      smooth(ring, values, threshold);
      if (area(ring) > 0) polygons.push([ring]);
      else holes.push(ring);
    });

    // Restore original dimensions
    dx = originalDx; dy = originalDy;

    // Enhanced hole assignment algorithm with spatial indexing
    holes = assignHolesOptimized(polygons, holes);

    // Adjust coordinates if boundaries were extended
    if (extendBoundaries) {
      polygons = adjustExtendedCoordinates(polygons);
    }

    const result = {
      type: "MultiPolygon",
      value: threshold,
      coordinates: polygons
    };

    if (exportGeoJSON) {
      result.properties = {
        value: threshold,
        level: level
      };
    }

    return result;
  }

  // Enhanced Marching Squares for level-based contours
  function isoringsForLevel(values, levelMap, threshold, targetLevel, callback) {
    var fragmentByStart = new Map(),
        fragmentByEnd = new Map(),
        x, y, t0, t1, t2, t3;

    // Special case for the first row (y = -1, t2 = t3 = 0).
    x = y = -1;
    t1 = isLevelAbove(levelMap[0], targetLevel);
    processCase(t1 << 1, x, y, [], stitch);
    while (++x < dx - 1) {
      t0 = t1, t1 = isLevelAbove(levelMap[x + 1], targetLevel);
      processCase(t0 | t1 << 1, x, y, [], stitch);
    }
    processCase(t1 << 0, x, y, [], stitch);

    // General case for the intermediate rows.
    while (++y < dy - 1) {
      x = -1;
      t1 = isLevelAbove(levelMap[y * dx + dx], targetLevel);
      t2 = isLevelAbove(levelMap[y * dx], targetLevel);
      processCase(t1 << 1 | t2 << 2, x, y, [], stitch);
      while (++x < dx - 1) {
        t0 = t1, t1 = isLevelAbove(levelMap[y * dx + dx + x + 1], targetLevel);
        t3 = t2, t2 = isLevelAbove(levelMap[y * dx + x + 1], targetLevel);
        
        const caseIndex = t0 | t1 << 1 | t2 << 2 | t3 << 3;
        const cellValues = [
          values[y * dx + x],
          values[y * dx + x + 1],
          values[y * dx + dx + x],
          values[y * dx + dx + x + 1]
        ];
        
        processCase(caseIndex, x, y, cellValues, stitch);
      }
      processCase(t1 | t2 << 3, x, y, [], stitch);
    }

    // Special case for the last row (y = dy - 1, t0 = t1 = 0).
    x = -1;
    t2 = isLevelAbove(levelMap[y * dx], targetLevel);
    processCase(t2 << 2, x, y, [], stitch);
    while (++x < dx - 1) {
      t3 = t2, t2 = isLevelAbove(levelMap[y * dx + x + 1], targetLevel);
      processCase(t2 << 2 | t3 << 3, x, y, [], stitch);
    }
    processCase(t2 << 3, x, y, [], stitch);

    function processCase(caseIndex, x, y, cellValues, stitch) {
      // Handle saddle point disambiguation for cases 5 and 10
      if (enableSaddleDisambiguation && (caseIndex === 5 || caseIndex === 10) && cellValues.length === 4) {
        const saddleCase = disambiguateSaddle(cellValues, threshold, caseIndex);
        saddleCase.forEach(stitch);
      } else {
        cases[caseIndex].forEach(stitch);
      }
    }

    function stitch(line) {
      var start = [line[0][0] + x, line[0][1] + y],
          end = [line[1][0] + x, line[1][1] + y],
          startIndex = index(start),
          endIndex = index(end),
          f, g;
      
      if (f = fragmentByEnd.get(startIndex)) {
        if (g = fragmentByStart.get(endIndex)) {
          fragmentByEnd.delete(f.end);
          fragmentByStart.delete(g.start);
          if (f === g) {
            f.ring.push(end);
            callback(f.ring);
          } else {
            const merged = {start: f.start, end: g.end, ring: f.ring.concat(g.ring)};
            fragmentByStart.set(f.start, merged);
            fragmentByEnd.set(g.end, merged);
          }
        } else {
          fragmentByEnd.delete(f.end);
          f.ring.push(end);
          fragmentByEnd.set(f.end = endIndex, f);
        }
      } else if (f = fragmentByStart.get(endIndex)) {
        if (g = fragmentByEnd.get(startIndex)) {
          fragmentByStart.delete(f.start);
          fragmentByEnd.delete(g.end);
          if (f === g) {
            f.ring.push(end);
            callback(f.ring);
          } else {
            const merged = {start: g.start, end: f.end, ring: g.ring.concat(f.ring)};
            fragmentByStart.set(g.start, merged);
            fragmentByEnd.set(f.end, merged);
          }
        } else {
          fragmentByStart.delete(f.start);
          f.ring.unshift(start);
          fragmentByStart.set(f.start = startIndex, f);
        }
      } else {
        const fragment = {start: startIndex, end: endIndex, ring: [start, end]};
        fragmentByStart.set(startIndex, fragment);
        fragmentByEnd.set(endIndex, fragment);
      }
    }
  }

  // Check if a level is above the target (for non-overlapping contours)
  function isLevelAbove(cellLevel, targetLevel) {
    if (cellLevel === -1) return false; // Invalid cell
    return cellLevel >= targetLevel;
  }

  // Extend grid boundaries to avoid edge artifacts
  function extendGridBoundaries(values) {
    const extended = new Array((dx + 2) * (dy + 2));
    
    // Fill extended grid
    for (let j = 0; j < dy + 2; j++) {
      for (let i = 0; i < dx + 2; i++) {
        const extIndex = j * (dx + 2) + i;
        
        if (i === 0 || i === dx + 1 || j === 0 || j === dy + 1) {
          // Boundary cells - extrapolate from nearest valid cell
          const nearestI = Math.max(0, Math.min(dx - 1, i - 1));
          const nearestJ = Math.max(0, Math.min(dy - 1, j - 1));
          const nearestIndex = nearestJ * dx + nearestI;
          extended[extIndex] = values[nearestIndex];
        } else {
          // Interior cells - copy from original
          const origIndex = (j - 1) * dx + (i - 1);
          extended[extIndex] = values[origIndex];
        }
      }
    }
    
    return extended;
  }

  // Adjust coordinates after boundary extension
  function adjustExtendedCoordinates(polygons) {
    return polygons.map(polygon => 
      polygon.map(ring => 
        ring.map(point => [point[0] - 1, point[1] - 1])
      )
    );
  }

  // Enhanced contour generation (fallback for single contours)
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

    // Enhanced hole assignment algorithm with spatial indexing
    holes = assignHolesOptimized(polygons, holes);

    const result = {
      type: "MultiPolygon",
      value: value,
      coordinates: polygons
    };

    if (exportGeoJSON) {
      result.properties = {
        value: value,
        level: value
      };
    }

    return result;
  }

  // Original isorings for single contour (kept for compatibility)
  function isorings(values, value, callback) {
    var fragmentByStart = new Map(),
        fragmentByEnd = new Map(),
        x, y, t0, t1, t2, t3;

    // Special case for the first row (y = -1, t2 = t3 = 0).
    x = y = -1;
    t1 = above(values[0], value);
    processCase(t1 << 1, x, y, [], stitch);
    while (++x < dx - 1) {
      t0 = t1, t1 = above(values[x + 1], value);
      processCase(t0 | t1 << 1, x, y, [], stitch);
    }
    processCase(t1 << 0, x, y, [], stitch);

    // General case for the intermediate rows.
    while (++y < dy - 1) {
      x = -1;
      t1 = above(values[y * dx + dx], value);
      t2 = above(values[y * dx], value);
      processCase(t1 << 1 | t2 << 2, x, y, [], stitch);
      while (++x < dx - 1) {
        t0 = t1, t1 = above(values[y * dx + dx + x + 1], value);
        t3 = t2, t2 = above(values[y * dx + x + 1], value);
        
        const caseIndex = t0 | t1 << 1 | t2 << 2 | t3 << 3;
        const cellValues = [
          values[y * dx + x],
          values[y * dx + x + 1],
          values[y * dx + dx + x],
          values[y * dx + dx + x + 1]
        ];
        
        processCase(caseIndex, x, y, cellValues, stitch);
      }
      processCase(t1 | t2 << 3, x, y, [], stitch);
    }

    // Special case for the last row (y = dy - 1, t0 = t1 = 0).
    x = -1;
    t2 = above(values[y * dx], value);
    processCase(t2 << 2, x, y, [], stitch);
    while (++x < dx - 1) {
      t3 = t2, t2 = above(values[y * dx + x + 1], value);
      processCase(t2 << 2 | t3 << 3, x, y, [], stitch);
    }
    processCase(t2 << 3, x, y, [], stitch);

    function processCase(caseIndex, x, y, cellValues, stitch) {
      // Handle saddle point disambiguation for cases 5 and 10
      if (enableSaddleDisambiguation && (caseIndex === 5 || caseIndex === 10) && cellValues.length === 4) {
        const saddleCase = disambiguateSaddle(cellValues, value, caseIndex);
        saddleCase.forEach(stitch);
      } else {
        cases[caseIndex].forEach(stitch);
      }
    }

    function stitch(line) {
      var start = [line[0][0] + x, line[0][1] + y],
          end = [line[1][0] + x, line[1][1] + y],
          startIndex = index(start),
          endIndex = index(end),
          f, g;
      
      if (f = fragmentByEnd.get(startIndex)) {
        if (g = fragmentByStart.get(endIndex)) {
          fragmentByEnd.delete(f.end);
          fragmentByStart.delete(g.start);
          if (f === g) {
            f.ring.push(end);
            callback(f.ring);
          } else {
            const merged = {start: f.start, end: g.end, ring: f.ring.concat(g.ring)};
            fragmentByStart.set(f.start, merged);
            fragmentByEnd.set(g.end, merged);
          }
        } else {
          fragmentByEnd.delete(f.end);
          f.ring.push(end);
          fragmentByEnd.set(f.end = endIndex, f);
        }
      } else if (f = fragmentByStart.get(endIndex)) {
        if (g = fragmentByEnd.get(startIndex)) {
          fragmentByStart.delete(f.start);
          fragmentByEnd.delete(g.end);
          if (f === g) {
            f.ring.push(end);
            callback(f.ring);
          } else {
            const merged = {start: g.start, end: f.end, ring: g.ring.concat(f.ring)};
            fragmentByStart.set(g.start, merged);
            fragmentByEnd.set(f.end, merged);
          }
        } else {
          fragmentByStart.delete(f.start);
          f.ring.unshift(start);
          fragmentByStart.set(f.start = startIndex, f);
        }
      } else {
        const fragment = {start: startIndex, end: endIndex, ring: [start, end]};
        fragmentByStart.set(startIndex, fragment);
        fragmentByEnd.set(endIndex, fragment);
      }
    }
  }

  // Saddle point disambiguation based on cell center value
  function disambiguateSaddle(cellValues, threshold, caseIndex) {
    const [v00, v01, v10, v11] = cellValues.map(v => validValue(v));
    const centerValue = (v00 + v01 + v10 + v11) / 4;
    
    // Choose case based on whether center is above or below threshold
    const useConnected = (centerValue >= threshold);
    const selectedCases = saddleCases[caseIndex];
    
    return selectedCases[useConnected ? 1 : 0];
  }

  // Optimized hole assignment with spatial indexing
  function assignHolesOptimized(polygons, holes) {
    if (holes.length === 0) return holes;
    
    // Create bounding boxes for polygons
    const polygonBounds = polygons.map(polygon => getBounds(polygon[0]));
    
    holes.forEach(function(hole) {
      const holeBounds = getBounds(hole);
      let assigned = false;
      
      // Find potential containing polygons using bounding box test
      for (let i = 0; i < polygons.length && !assigned; i++) {
        if (boundsContain(polygonBounds[i], holeBounds)) {
          if (contains(polygons[i][0], hole) !== -1) {
            polygons[i].push(hole);
            assigned = true;
          }
        }
      }
    });
    
    return [];
  }

  // Calculate bounding box for a ring
  function getBounds(ring) {
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (const [x, y] of ring) {
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    }
    return {minX, minY, maxX, maxY};
  }

  // Check if bounds1 contains bounds2
  function boundsContain(bounds1, bounds2) {
    return bounds1.minX <= bounds2.minX && bounds1.maxX >= bounds2.maxX &&
           bounds1.minY <= bounds2.minY && bounds1.maxY >= bounds2.maxY;
  }

  function index(point) {
    return point[0] * 2 + point[1] * (dx + 1) * 4;
  }

  // Enhanced smoothing with adjustable factor
  function smoothLinear(ring, values, value) {
    ring.forEach(function(point) {
      var x = point[0],
          y = point[1],
          xt = x | 0,
          yt = y | 0,
          v1 = validValue(values[yt * dx + xt]);
      
      if (x > 0 && x < dx && xt === x) {
        const smoothed = smooth1(x, validValue(values[yt * dx + xt - 1]), v1, value);
        point[0] = x + (smoothed - x) * smoothFactor;
      }
      if (y > 0 && y < dy && yt === y) {
        const smoothed = smooth1(y, validValue(values[(yt - 1) * dx + xt]), v1, value);
        point[1] = y + (smoothed - y) * smoothFactor;
      }
    });
  }

  // Enhanced API
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

  contours.smoothFactor = function(_) {
    return arguments.length ? (smoothFactor = Math.max(0, Math.min(1, +_)), contours) : smoothFactor;
  };

  contours.saddleDisambiguation = function(_) {
    return arguments.length ? (enableSaddleDisambiguation = !!_, contours) : enableSaddleDisambiguation;
  };

  contours.geoJSON = function(_) {
    return arguments.length ? (exportGeoJSON = !!_, contours) : exportGeoJSON;
  };

  contours.nullHandling = function(_) {
    return arguments.length ? (handleNullValues = !!_, contours) : handleNullValues;
  };

  // New API for topology control
  contours.preventOverlap = function(_) {
    return arguments.length ? (preventOverlap = !!_, contours) : preventOverlap;
  };

  contours.extendBoundaries = function(_) {
    return arguments.length ? (extendBoundaries = !!_, contours) : extendBoundaries;
  };

  return contours;
}

// When computing the extent, ignore infinite values (as well as invalid ones).
function finite(x) {
  return isFinite(x) ? x : NaN;
}

// Enhanced null value handling
function above(x, value) {
  if (x == null) return false;
  const num = +x;
  return isFinite(num) ? num >= value : false;
}

// Enhanced value validation with null handling
function validValue(v) {
  if (v == null) return -Infinity;
  const num = +v;
  return isFinite(num) ? num : -Infinity;
}

function smooth1(x, v0, v1, value) {
  const a = value - v0;
  const b = v1 - v0;
  const d = isFinite(a) && isFinite(b) && b !== 0 ? a / b : 0.5;
  return isNaN(d) ? x : x + d - 0.5;
}