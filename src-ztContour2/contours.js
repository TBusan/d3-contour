import {extent, nice, thresholdSturges, ticks} from "d3-array";
import {
  area, 
  ascending, 
  above, 
  constant, 
  createSmoothingFunction,
  hasSelfIntersections, 
  SpatialIndex
} from "./utils.js";

// Define the marching squares cases as in the original d3-contour
// Each case is an array of line segments connecting edges of a unit square
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

export default function() {
  var dx = 1,
      dy = 1,
      threshold = thresholdSturges,
      smoothing = 1.0,
      mode = "surfaces",  // "surfaces", "lines", or "both"
      nullValue = null,   // Value to treat as null/missing
      chunkSize = 10000,  // Process data in chunks to reduce memory pressure
      skipHoles = false;  // Skip hole processing for complex datasets with many null values

  function contours(values) {
    // Handle null values by replacing them with NaN
    if (nullValue !== undefined) {
      values = handleNullValues(values);
    }
    
    var tz = computeThresholds(values);
    
    // Limit number of thresholds for large datasets
    if (values.length > 1000000 && tz.length > 10) {
      console.warn("Large dataset detected, limiting number of thresholds to prevent memory issues");
      tz = tz.slice(0, 10);
    }

    // Generate contours based on the selected mode
    if (mode === "lines") {
      return tz.map(value => contourLines(values, value));
    } else if (mode === "both") {
      return tz.map(value => ({
        value: value,
        surfaces: contourSurfaces(values, value),
        lines: contourLines(values, value)
      }));
    } else { // Default: surfaces
      return tz.map(value => contourSurfaces(values, value));
    }
  }
  
  // Process null values in the input data
  function handleNullValues(values) {
    const n = values.length;
    
    // Check if we're dealing with a sparse dataset with many nulls
    let nullCount = 0;
    for (let i = 0; i < Math.min(n, 1000); i++) {
      const val = values[i];
      if (val === nullValue || val === null || val === undefined) {
        nullCount++;
      }
    }
    
    // If more than 30% of the sampled values are null, use the skipHoles option
    if (nullCount > 300) {
      skipHoles = true;
    }
    
    // For very large arrays, process in chunks to avoid memory pressure
    if (n > chunkSize) {
      const result = new Float64Array(n);
      
      for (let i = 0; i < n; i += chunkSize) {
        const chunk = Math.min(chunkSize, n - i);
        for (let j = 0; j < chunk; j++) {
          const val = values[i + j];
          result[i + j] = val === nullValue || val === null || val === undefined ? NaN : val;
        }
      }
      
      return result;
    } else {
      const result = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        const val = values[i];
        result[i] = val === nullValue || val === null || val === undefined ? NaN : val;
      }
      return result;
    }
  }

  // Compute threshold values for contours
  function computeThresholds(values) {
    let tz = threshold(values);
    
    // Convert number of thresholds into uniform thresholds
    if (!Array.isArray(tz)) {
      // Use sampling for large datasets to calculate extent
      let e;
      if (values.length > 100000) {
        const sampledValues = [];
        const sampleRate = Math.max(1, Math.floor(values.length / 10000));
        for (let i = 0; i < values.length; i += sampleRate) {
          const val = values[i];
          if (val !== null && val !== undefined && !isNaN(val) && isFinite(val)) {
            sampledValues.push(val);
          }
        }
        e = extent(sampledValues);
      } else {
        // Filter out null/NaN values before computing extent
        const validValues = [];
        for (let i = 0; i < values.length; i++) {
          const val = values[i];
          if (val !== null && val !== undefined && !isNaN(val) && isFinite(val)) {
            validValues.push(val);
          }
        }
        e = extent(validValues);
      }
      
      if (e[0] === undefined || e[1] === undefined || e[0] === e[1]) {
        console.warn("Invalid data range detected, returning empty contours");
        return [];
      }
      
      try {
        tz = ticks(...nice(e[0], e[1], tz), tz);
        
        // Filter out thresholds outside data range
        while (tz.length > 0 && tz[tz.length - 1] >= e[1]) tz.pop();
        while (tz.length > 0 && tz[0] < e[0]) tz.shift();
      } catch (err) {
        console.warn("Error calculating thresholds:", err.message);
        // Fallback: create simple thresholds
        const range = e[1] - e[0];
        tz = [];
        for (let i = 1; i < 10; i++) {
          tz.push(e[0] + (range * i / 10));
        }
      }
    } else {
      tz = tz.slice().sort(ascending);
    }
    
    return tz;
  }

  // Enhanced contour generator for surfaces (MultiPolygon)
  function contourSurfaces(values, value) {
    const v = value == null ? NaN : +value;
    if (isNaN(v)) throw new Error(`invalid value: ${value}`);

    var polygons = [],
        holes = [];

    // Generate contour rings
    isorings(values, v, function(ring) {
      // Skip tiny rings that are likely artifacts from noisy data
      if (ring.length < 4) return;
      
      // Apply smoothing if enabled
      if (smoothing > 0) {
        const smoothFunc = createSmoothingFunction(smoothing);
        smoothFunc(ring, values, v, dx, dy);
      }
      
      // Check if ring is a polygon (clockwise) or hole (counter-clockwise)
      const ringArea = area(ring);
      
      // Skip rings with zero area
      if (Math.abs(ringArea) < 1e-6) return;
      
      if (ringArea > 0) {
        // For large datasets, skip self-intersection checks to save memory
        if (values.length < 250000 && hasSelfIntersections(ring)) {
          console.warn("Self-intersecting ring detected");
        }
        polygons.push([ring]);
      } else {
        // Only collect holes if we're not skipping hole processing
        if (!skipHoles) {
          holes.push(ring);
        }
      }
    });

    // For datasets with many nulls, we might want to skip hole assignment
    // as it can lead to incorrect topology
    if (!skipHoles && holes.length > 0 && polygons.length > 0) {
      if (holes.length > 1000 || polygons.length > 1000) {
        // For extremely large datasets, use a simplified approach
        assignHolesSimplified(holes, polygons);
      } else {
        const spatialIndex = new SpatialIndex();
        
        // Add all polygons to the spatial index
        polygons.forEach((polygon, i) => {
          spatialIndex.addPolygon(polygon, i);
        });
        
        // Assign each hole to the smallest containing polygon
        holes.forEach(function(hole) {
          const idx = spatialIndex.findContainingPolygon(hole);
          if (idx !== -1) {
            polygons[idx].push(hole);
          }
        });
      }
    }

    return {
      type: "MultiPolygon",
      value: value,
      coordinates: polygons
    };
  }
  
  // Simplified hole assignment for very large datasets
  function assignHolesSimplified(holes, polygons) {
    // Use centroids and bounding boxes for quick assignment
    const polyBounds = polygons.map(polygon => {
      const ring = polygon[0];
      let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
      
      for (const [x, y] of ring) {
        if (x < minX) minX = x;
        if (y < minY) minY = y;
        if (x > maxX) maxX = x;
        if (y > maxY) maxY = y;
      }
      
      return {minX, minY, maxX, maxY};
    });
    
    // Process holes in batches
    const batchSize = 200;
    for (let i = 0; i < holes.length; i += batchSize) {
      const batch = holes.slice(i, i + batchSize);
      
      batch.forEach(hole => {
        // Find the centroid of the hole
        let cx = 0, cy = 0;
        for (const [x, y] of hole) {
          cx += x;
          cy += y;
        }
        cx /= hole.length;
        cy /= hole.length;
        
        // Find containing polygon using bounding box test
        for (let j = 0; j < polygons.length; j++) {
          const bounds = polyBounds[j];
          if (cx >= bounds.minX && cx <= bounds.maxX && cy >= bounds.minY && cy <= bounds.maxY) {
            // Simplified point-in-polygon test
            const ring = polygons[j][0];
            let inside = false;
            
            for (let k = 0, l = ring.length - 1; k < ring.length; l = k++) {
              const [xi, yi] = ring[k];
              const [xj, yj] = ring[l];
              
              if (((yi > cy) !== (yj > cy)) && 
                  (cx < (xj - xi) * (cy - yi) / (yj - yi) + xi)) {
                inside = !inside;
              }
            }
            
            if (inside) {
              polygons[j].push(hole);
              break;
            }
          }
        }
      });
    }
  }
  
  // Enhanced contour generator for lines (MultiLineString)
  function contourLines(values, value) {
    const v = value == null ? NaN : +value;
    if (isNaN(v)) throw new Error(`invalid value: ${value}`);
    
    const lines = [];
    
    // Generate contour lines
    isoRingsAsLines(values, v, function(line) {
      // Skip tiny lines that are likely artifacts
      if (line.length < 3) return;
      
      // Apply smoothing if enabled (skip for very large datasets)
      if (smoothing > 0 && values.length < 250000) {
        const smoothFunc = createSmoothingFunction(smoothing);
        smoothFunc(line, values, v, dx, dy);
      }
      
      lines.push(line);
    });
    
    return {
      type: "MultiLineString",
      value: value,
      coordinates: lines
    };
  }

  // Generate iso-value rings (polygons)
  // Memory-optimized version of the original algorithm
  function isorings(values, value, callback) {
    var fragmentByStart = new Map(),
        fragmentByEnd = new Map(),
        fragments = [];  // Keep track of all fragments for post-processing

    // Trace through the grid using marching squares
    marchingSquares(values, value, function(x, y, caseIndex) {
      // Skip invalid cases (for sparse data with many nulls)
      if (isNaN(x) || isNaN(y)) return;
      
      const segments = cases[caseIndex];
      for (const segment of segments) {
        const start = [segment[0][0] + x, segment[0][1] + y];
        const end = [segment[1][0] + x, segment[1][1] + y];
        
        const startKey = index(start);
        const endKey = index(end);
        
        let f = fragmentByEnd.get(startKey);
        let g = fragmentByStart.get(endKey);
        
        // Various fragment connection cases
        if (f !== undefined && g !== undefined) {
          fragmentByEnd.delete(startKey);
          fragmentByStart.delete(endKey);
          
          if (f === g) {
            // We've closed a loop
            f.ring.push(end);
            
            // Check for minimum length and self-intersections
            if (f.ring.length >= 4) {
              // Remove self-intersections before returning the ring
              const cleanedRing = removeSelfIntersections(f.ring);
              callback(cleanedRing);
            }
            
            // Clean up references to allow garbage collection
            f.ring = null;
          } else {
            // Connect two fragments
            // Use object pooling to reduce allocations
            const combinedRing = f.ring.concat(g.ring);
            f.ring = null;
            g.ring = null;
            
            const newFragment = {
              startKey: f.startKey,
              endKey: g.endKey,
              ring: combinedRing
            };
            
            fragmentByStart.set(f.startKey, newFragment);
            fragmentByEnd.set(g.endKey, newFragment);
            
            // Replace in fragments list
            const fIndex = fragments.indexOf(f);
            const gIndex = fragments.indexOf(g);
            if (fIndex !== -1) fragments[fIndex] = newFragment;
            if (gIndex !== -1) fragments.splice(gIndex, 1);
            else fragments.push(newFragment);
          }
        } else if (f !== undefined) {
          // Extend existing fragment from end
          fragmentByEnd.delete(startKey);
          f.ring.push(end);
          f.endKey = endKey;
          fragmentByEnd.set(endKey, f);
        } else if (g !== undefined) {
          // Extend existing fragment from start
          fragmentByStart.delete(endKey);
          g.ring.unshift(start);
          g.startKey = startKey;
          fragmentByStart.set(startKey, g);
        } else {
          // Create a new fragment
          const fragment = {
            startKey: startKey,
            endKey: endKey,
            ring: [start, end]
          };
          fragmentByStart.set(startKey, fragment);
          fragmentByEnd.set(endKey, fragment);
          fragments.push(fragment);
        }
      }
    });
    
    // Attempt to join nearby open fragments first
    attemptToJoinFragments(fragments, fragmentByStart, fragmentByEnd);
    
    // Process any remaining fragments (these would be open contours at the boundary)
    // For very large maps, skip boundary closing to save memory
    if (values.length < 250000) {
      fragmentByStart.forEach(function(fragment, key) {
        // Skip fragments that are already closed
        if (fragmentByEnd.has(key)) return;
        
        // Close open contours at boundaries if needed
        if (isBoundaryFragment(fragment.ring, dx, dy)) {
          closeBoundaryFragment(fragment.ring, dx, dy);
          
          // Remove any self-intersections that might have been created by closing
          const cleanedRing = removeSelfIntersections(fragment.ring);
          if (cleanedRing.length >= 4) {
            callback(cleanedRing);
          }
        } else if (fragment.ring.length >= 4) {
          // Try to close the ring by connecting start and end if they're close
          const start = fragment.ring[0];
          const end = fragment.ring[fragment.ring.length - 1];
          
          if (getDistance(start, end) < Math.sqrt(dx*dx + dy*dy) / 2) {
            fragment.ring.push([...start]); // Close the loop by adding a copy of the start
            callback(fragment.ring);
          } else {
            // For lines within the grid that couldn't be closed,
            // attempt to interpolate between sparse data regions
            interpolateOpenFragment(fragment.ring, values, value);
            callback(fragment.ring);
          }
        }
      });
    } else {
      // For large datasets with sparse data, try to close fragments if possible
      fragmentByStart.forEach(function(fragment, key) {
        // Skip fragments that are already closed
        if (fragmentByEnd.has(key)) return;
        
        // If the fragment has enough points, consider it valid
        if (fragment.ring.length >= 4) {
          // Simple closing - just connect the endpoints if they're close
          const start = fragment.ring[0];
          const end = fragment.ring[fragment.ring.length - 1];
          
          if (getDistance(start, end) < Math.sqrt(dx*dx + dy*dy)) {
            fragment.ring.push([...start]);
          }
          callback(fragment.ring);
        }
      });
    }
    
    // Clear maps to help garbage collection
    fragmentByStart.clear();
    fragmentByEnd.clear();
    fragments.length = 0;
  }
  
  // Attempt to join nearby open fragments
  function attemptToJoinFragments(fragments, fragmentByStart, fragmentByEnd) {
    // For each open fragment, check if another fragment's start/end point is close to its end/start
    let joined;
    do {
      joined = false;
      
      for (let i = 0; i < fragments.length; i++) {
        const f = fragments[i];
        
        // Skip if the fragment is no longer open
        if (!fragmentByStart.has(f.startKey) || !fragmentByEnd.has(f.endKey)) continue;
        
        const fEnd = f.ring[f.ring.length - 1];
        
        // Check all other fragments
        for (let j = 0; j < fragments.length; j++) {
          if (i === j) continue;
          
          const g = fragments[j];
          
          // Skip if the fragment is no longer open
          if (!fragmentByStart.has(g.startKey) || !fragmentByEnd.has(g.endKey)) continue;
          
          const gStart = g.ring[0];
          
          // If endpoints are close, join the fragments
          if (getDistance(fEnd, gStart) < 0.5) {
            // Remove the end point from maps
            fragmentByEnd.delete(f.endKey);
            fragmentByStart.delete(g.startKey);
            
            // Join the rings
            const combinedRing = f.ring.concat(g.ring);
            
            // Create new fragment
            const newFragment = {
              startKey: f.startKey,
              endKey: g.endKey,
              ring: combinedRing
            };
            
            // Update maps
            fragmentByStart.set(f.startKey, newFragment);
            fragmentByEnd.set(g.endKey, newFragment);
            
            // Replace in fragments list
            fragments[i] = newFragment;
            fragments.splice(j, 1);
            
            joined = true;
            break;
          }
        }
        
        if (joined) break;
      }
    } while (joined);
  }
  
  // Interpolate points for open fragments that might be caused by sparse data regions
  function interpolateOpenFragment(ring, values, value) {
    // Check if the fragment has endpoints in valid data regions
    const start = ring[0];
    const end = ring[ring.length - 1];
    
    // If the endpoints are too far away, don't try to interpolate
    if (getDistance(start, end) > Math.sqrt(dx*dx + dy*dy) * 3) {
      return;
    }
    
    // Create a straight line between start and end
    const steps = 10; // Number of interpolation steps
    const dx = (end[0] - start[0]) / steps;
    const dy = (end[1] - start[1]) / steps;
    
    // Add interpolated points
    for (let i = 1; i < steps; i++) {
      ring.push([
        start[0] + dx * i,
        start[1] + dy * i
      ]);
    }
    
    // Close the loop by adding the first point again
    ring.push([...ring[0]]);
  }
  
  // Remove self-intersections from a ring by simplifying problematic segments
  function removeSelfIntersections(ring) {
    if (ring.length < 4) return ring;
    
    const result = [...ring];
    
    // Check for self-intersections
    if (!hasSelfIntersections(result)) {
      return result;
    }
    
    // Find and resolve self-intersections
    const newRing = [];
    const used = new Array(result.length).fill(false);
    
    // Start with the first point
    newRing.push(result[0]);
    used[0] = true;
    
    let current = 0;
    
    while (newRing.length < result.length) {
      // Find the closest unused point that doesn't create an intersection
      let bestDist = Infinity;
      let bestIdx = -1;
      
      for (let i = 0; i < result.length; i++) {
        if (used[i]) continue;
        
        // Check if adding this point would create a self-intersection
        const wouldIntersect = wouldCreateIntersection(newRing, result[i]);
        
        if (!wouldIntersect) {
          const dist = getDistance(result[current], result[i]);
          if (dist < bestDist) {
            bestDist = dist;
            bestIdx = i;
          }
        }
      }
      
      // If we found a valid next point, add it
      if (bestIdx !== -1) {
        newRing.push(result[bestIdx]);
        used[bestIdx] = true;
        current = bestIdx;
      } else {
        // If we can't find a valid next point, try to close the ring
        if (!used[0] && !wouldCreateIntersection(newRing, result[0])) {
          newRing.push(result[0]);
        }
        break;
      }
    }
    
    // If the resulting ring is too small, return the original
    if (newRing.length < 4) {
      return ring;
    }
    
    return newRing;
  }
  
  // Check if adding a point to a path would create a self-intersection
  function wouldCreateIntersection(path, point) {
    if (path.length < 3) return false;
    
    const last = path[path.length - 1];
    const newSegment = [last, point];
    
    // Check against all non-adjacent segments
    for (let i = 0; i < path.length - 2; i++) {
      const segment = [path[i], path[i + 1]];
      if (segmentsIntersect(newSegment[0], newSegment[1], segment[0], segment[1])) {
        return true;
      }
    }
    
    return false;
  }
  
  // Check if two line segments intersect
  function segmentsIntersect(a, b, c, d) {
    // Check if segments share an endpoint
    if (pointsEqual(a, c) || pointsEqual(a, d) || pointsEqual(b, c) || pointsEqual(b, d)) {
      return false;
    }
    
    const ccw = (p1, p2, p3) => {
      return (p3[1] - p1[1]) * (p2[0] - p1[0]) > (p2[1] - p1[1]) * (p3[0] - p1[0]);
    };
    
    return ccw(a, c, d) !== ccw(b, c, d) && ccw(a, b, c) !== ccw(a, b, d);
  }
  
  function pointsEqual(p1, p2) {
    return p1[0] === p2[0] && p1[1] === p2[1];
  }

  // Memory-optimized version for line generation
  function isoRingsAsLines(values, value, callback) {
    var fragmentByStart = new Map(),
        fragmentByEnd = new Map();

    // Trace through the grid using marching squares
    marchingSquares(values, value, function(x, y, caseIndex) {
      // Skip invalid cases (for sparse data with many nulls)
      if (isNaN(x) || isNaN(y)) return;
      
      const segments = cases[caseIndex];
      for (const segment of segments) {
        const start = [segment[0][0] + x, segment[0][1] + y];
        const end = [segment[1][0] + x, segment[1][1] + y];
        
        const startKey = index(start);
        const endKey = index(end);
        
        let f = fragmentByEnd.get(startKey);
        let g = fragmentByStart.get(endKey);
        
        // Various fragment connection cases
        if (f !== undefined && g !== undefined) {
          fragmentByEnd.delete(startKey);
          fragmentByStart.delete(endKey);
          
          if (f === g) {
            // We've closed a loop - for lines we still want to preserve the loop
            f.ring.push(end);
            callback(f.ring);
            f.ring = null; // Help garbage collection
          } else {
            // Connect two fragments
            const connectedRing = f.ring.concat(g.ring);
            f.ring = null;
            g.ring = null;
            
            fragmentByStart.set(f.startKey, {
              startKey: f.startKey,
              endKey: g.endKey,
              ring: connectedRing
            });
            fragmentByEnd.set(g.endKey, fragmentByStart.get(f.startKey));
          }
        } else if (f !== undefined) {
          // Extend existing fragment from end
          fragmentByEnd.delete(startKey);
          f.ring.push(end);
          f.endKey = endKey;
          fragmentByEnd.set(endKey, f);
        } else if (g !== undefined) {
          // Extend existing fragment from start
          fragmentByStart.delete(endKey);
          g.ring.unshift(start);
          g.startKey = startKey;
          fragmentByStart.set(startKey, g);
        } else {
          // Create a new fragment
          const fragment = {
            startKey: startKey,
            endKey: endKey,
            ring: [start, end]
          };
          fragmentByStart.set(startKey, fragment);
          fragmentByEnd.set(endKey, fragment);
        }
      }
    });
    
    // For lines, we want to return all fragments with adequate length
    fragmentByStart.forEach(function(fragment, key) {
      if (!fragmentByEnd.has(key) && fragment.ring.length >= 3) {
        callback(fragment.ring);
        fragment.ring = null; // Help garbage collection
      }
    });
    
    // Clear maps to help garbage collection
    fragmentByStart.clear();
    fragmentByEnd.clear();
  }
  
  // Check if a fragment touches the boundary
  function isBoundaryFragment(ring, dx, dy) {
    for (const [x, y] of ring) {
      if (x <= 0 || x >= dx || y <= 0 || y >= dy) {
        return true;
      }
    }
    return false;
  }
  
  // Close an open boundary fragment
  function closeBoundaryFragment(ring, dx, dy) {
    const start = ring[0];
    const end = ring[ring.length - 1];
    
    // If already closed, do nothing
    if (start[0] === end[0] && start[1] === end[1]) return;
    
    // Find closest boundary points - use simple approach to save memory
    let startBoundary = [
      Math.max(0, Math.min(dx, start[0])),
      Math.max(0, Math.min(dy, start[1]))
    ];
    
    let endBoundary = [
      Math.max(0, Math.min(dx, end[0])),
      Math.max(0, Math.min(dy, end[1]))
    ];
    
    // Simplified boundary closing - just connect directly
    ring.push(endBoundary, startBoundary);
  }

  // Memory-optimized marching squares implementation
  // Enhanced to better handle sparse data with null/NaN values
  function marchingSquares(values, value, callback) {
    // Process in rows to reduce memory pressure
    const processRow = (y) => {
      if (y < -1 || y >= dy) return; // Skip out-of-bounds rows
      
      let x = -1;
      let t1, t2, t3;
      
      if (y === -1) {
        // Special case for the first row (y = -1)
        t1 = isValidForContour(values[0]) ? above(values[0], value) : false;
        cases[t1 << 1].forEach(() => callback(x, y, t1 << 1));
        
        while (++x < dx - 1) {
          t1 = isValidForContour(values[x + 1]) ? above(values[x + 1], value) : false;
          cases[t1 << 1].forEach(() => callback(x, y, t1 << 1));
        }
        
        cases[t1 << 0].forEach(() => callback(x, y, t1 << 0));
      } else if (y === dy - 1) {
        // Special case for the last row
        x = -1;
        t2 = isValidForContour(values[y * dx]) ? above(values[y * dx], value) : false;
        cases[t2 << 2].forEach(() => callback(x, y, t2 << 2));
        
        while (++x < dx - 1) {
          t3 = t2;
          t2 = isValidForContour(values[y * dx + x + 1]) ? above(values[y * dx + x + 1], value) : false;
          cases[t2 << 2 | t3 << 3].forEach(() => callback(x, y, t2 << 2 | t3 << 3));
        }
        
        cases[t2 << 3].forEach(() => callback(x, y, t2 << 3));
      } else {
        // General case for intermediate rows
        x = -1;
        t1 = isValidForContour(values[y * dx + dx]) ? above(values[y * dx + dx], value) : false;
        t2 = isValidForContour(values[y * dx]) ? above(values[y * dx], value) : false;
        cases[t1 << 1 | t2 << 2].forEach(() => callback(x, y, t1 << 1 | t2 << 2));
        
        while (++x < dx - 1) {
          const t0 = t1;
          t1 = isValidForContour(values[y * dx + dx + x + 1]) ? above(values[y * dx + dx + x + 1], value) : false;
          const t3 = t2;
          t2 = isValidForContour(values[y * dx + x + 1]) ? above(values[y * dx + x + 1], value) : false;
          
          // Skip cells with NaN or null values
          if (
            !isValidForContour(values[y * dx + x]) || 
            !isValidForContour(values[y * dx + x + 1]) || 
            !isValidForContour(values[y * dx + dx + x]) || 
            !isValidForContour(values[y * dx + dx + x + 1])
          ) {
            // If some corners have valid values, try to interpolate
            tryInterpolateCell(x, y, values, value, callback);
            continue;
          }
          
          let caseIndex = t0 | t1 << 1 | t2 << 2 | t3 << 3;
          
          // Always apply saddle point disambiguation for better topological correctness
          // For cases 5 and 10, use advanced disambiguation
          if (caseIndex === 5 || caseIndex === 10) {
            // Get corner values for saddle point disambiguation
            const corners = [
              values[y * dx + x],           // bottom-left
              values[y * dx + x + 1],       // bottom-right
              values[y * dx + dx + x + 1],  // top-right
              values[y * dx + dx + x]       // top-left
            ];
            
            // Use enhanced disambiguation with bilinear interpolation for better accuracy
            caseIndex = disambiguateSaddleEnhanced(caseIndex, corners, value, x, y);
          }
          
          // Create segments for this case
          cases[caseIndex].forEach(() => callback(x, y, caseIndex));
        }
        
        cases[t1 | t2 << 3].forEach(() => callback(x, y, t1 | t2 << 3));
      }
    };
    
    // Process rows in chunks for very large datasets
    const rowChunkSize = Math.min(100, dy);
    
    // Process first row (y = -1)
    processRow(-1);
    
    // Process intermediate rows in chunks
    for (let yChunk = 0; yChunk < dy - 1; yChunk += rowChunkSize) {
      const endY = Math.min(yChunk + rowChunkSize, dy - 1);
      for (let y = yChunk; y < endY; y++) {
        processRow(y);
      }
    }
    
    // Process last row if not already processed
    processRow(dy - 1);
  }
  
  // Try to interpolate contour line segments for cells with missing data
  function tryInterpolateCell(x, y, values, value, callback) {
    // Get the cell's corner values
    const i00 = y * dx + x;          // bottom-left
    const i10 = y * dx + x + 1;      // bottom-right
    const i11 = (y+1) * dx + x + 1;  // top-right
    const i01 = (y+1) * dx + x;      // top-left
    
    const v00 = isValidForContour(values[i00]) ? values[i00] : null;
    const v10 = isValidForContour(values[i10]) ? values[i10] : null;
    const v11 = isValidForContour(values[i11]) ? values[i11] : null;
    const v01 = isValidForContour(values[i01]) ? values[i01] : null;
    
    // Count valid corners
    let validCount = 0;
    if (v00 !== null) validCount++;
    if (v10 !== null) validCount++;
    if (v11 !== null) validCount++;
    if (v01 !== null) validCount++;
    
    // Not enough data to interpolate
    if (validCount < 2) return;
    
    // Check for value crossings and create interpolated line segments
    
    // Helper to interpolate position based on value
    function interpolateEdge(x0, y0, v0, x1, y1, v1, isoValue) {
      if (v0 === null || v1 === null) return null;
      if ((v0 >= isoValue && v1 >= isoValue) || (v0 < isoValue && v1 < isoValue)) return null;
      
      const t = (isoValue - v0) / (v1 - v0);
      return [x0 + t * (x1 - x0), y0 + t * (y1 - y0)];
    }
    
    // Check each edge
    // Bottom edge (v00-v10)
    const bottomPoint = interpolateEdge(x, y, v00, x+1, y, v10, value);
    
    // Right edge (v10-v11)
    const rightPoint = interpolateEdge(x+1, y, v10, x+1, y+1, v11, value);
    
    // Top edge (v11-v01)
    const topPoint = interpolateEdge(x+1, y+1, v11, x, y+1, v01, value);
    
    // Left edge (v01-v00)
    const leftPoint = interpolateEdge(x, y+1, v01, x, y, v00, value);
    
    // Collect valid points
    const points = [bottomPoint, rightPoint, topPoint, leftPoint].filter(p => p !== null);
    
    // If we have exactly 2 points, we can create a line segment
    if (points.length === 2) {
      // Create a custom case with a line segment between the two points
      const customCase = [[points[0], points[1]]];
      customCase.forEach(() => callback(x, y, -1)); // Use -1 to indicate custom case
    } else if (points.length > 2) {
      // For more than 2 points, create segments connecting adjacent points
      for (let i = 0; i < points.length; i++) {
        const j = (i + 1) % points.length;
        const customCase = [[points[i], points[j]]];
        customCase.forEach(() => callback(x, y, -1)); // Use -1 to indicate custom case
      }
    }
  }
  
  // Enhanced saddle point disambiguation using bilinear interpolation
  function disambiguateSaddleEnhanced(caseIndex, corners, value, x, y) {
    if (caseIndex !== 5 && caseIndex !== 10) return caseIndex;
    
    // For cases 5 and 10, use bilinear interpolation to determine the saddle configuration
    
    // First, handle traditional average-based approach as a fallback
    const avg = (corners[0] + corners[1] + corners[2] + corners[3]) / 4;
    
    // If the difference is significant, use the traditional approach
    if (Math.abs(avg - value) > Math.abs(corners[0] - corners[3]) * 0.1) {
      return value >= avg ? caseIndex : (caseIndex === 5 ? 10 : 5);
    }
    
    // Use bilinear interpolation to find more accurate center value
    // v(x,y) = a + b*x + c*y + d*x*y where (x,y) are in [0,1]
    const v00 = corners[0]; // bottom-left
    const v10 = corners[1]; // bottom-right
    const v11 = corners[2]; // top-right
    const v01 = corners[3]; // top-left
    
    // Bilinear interpolation coefficients
    const a = v00;
    const b = v10 - v00;
    const c = v01 - v00;
    const d = v11 - v01 - v10 + v00;
    
    // Find where the contour crosses the cell diagonals
    // Diagonal 1: (0,0) to (1,1)
    // t: parameter along diagonal where v(t,t) = value
    function solveQuadratic(a, b, c) {
      const discriminant = b * b - 4 * a * c;
      if (discriminant < 0) return [];
      if (discriminant === 0) return [-b / (2 * a)];
      const sqrtDisc = Math.sqrt(discriminant);
      return [(-b + sqrtDisc) / (2 * a), (-b - sqrtDisc) / (2 * a)];
    }
    
    // For diagonal 1: v(t,t) = a + b*t + c*t + d*t*t = value
    // This gives: a + (b+c)*t + d*t*t = value
    const diagonal1 = solveQuadratic(d, b + c, a - value);
    
    // For diagonal 2: v(t,1-t) = a + b*t + c*(1-t) + d*t*(1-t) = value
    // This gives: (a+c) + (b-c)*t + d*t*(1-t) = value
    // Expanded: (a+c) + (b-c)*t + d*t - d*t*t = value
    const diagonal2 = solveQuadratic(-d, b - c + d, a + c - value);
    
    // Filter solutions to those within cell [0,1]
    const valid1 = diagonal1.filter(t => t >= 0 && t <= 1);
    const valid2 = diagonal2.filter(t => t >= 0 && t <= 1);
    
    // If we found valid intersections, use them to determine case
    if (valid1.length > 0 || valid2.length > 0) {
      // If diagonal 1 has a valid intersection, use case 5
      if (valid1.length > 0) {
        return 5;
      }
      // If diagonal 2 has a valid intersection, use case 10
      if (valid2.length > 0) {
        return 10;
      }
    }
    
    // Fallback to traditional average-based approach
    return value >= avg ? caseIndex : (caseIndex === 5 ? 10 : 5);
  }

  // Helper function to check if a value is valid for contour generation
  function isValidForContour(val) {
    return val !== null && val !== undefined && !isNaN(val) && isFinite(val);
  }

  // Optimized helper function to create a unique index for a point
  function index(point) {
    return Math.floor(point[0] * 1000) + Math.floor(point[1] * 1000) * 2000000;
  }

  // Add a distance function
  function getDistance(p1, p2) {
    return Math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2);
  }

  // Public API
  contours.contour = contourSurfaces;
  
  contours.size = function(_) {
    if (!arguments.length) return [dx, dy];
    var _0 = Math.floor(_[0]), _1 = Math.floor(_[1]);
    if (!(_0 >= 0 && _1 >= 0)) throw new Error("invalid size");
    return dx = _0, dy = _1, contours;
  };
  
  contours.thresholds = function(_) {
    return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? constant(_.slice()) : constant(_), contours) : threshold;
  };
  
  contours.smooth = function(_) {
    return arguments.length ? (smoothing = _ === true ? 1 : _ === false ? 0 : +_, contours) : smoothing;
  };
  
  contours.mode = function(_) {
    return arguments.length ? (mode = _ + "", contours) : mode;
  };
  
  contours.nullValue = function(_) {
    return arguments.length ? (nullValue = _, contours) : nullValue;
  };
  
  contours.chunkSize = function(_) {
    return arguments.length ? (chunkSize = +_, contours) : chunkSize;
  };
  
  contours.skipHoles = function(_) {
    return arguments.length ? (skipHoles = !!_, contours) : skipHoles;
  };
  
  return contours;
} 