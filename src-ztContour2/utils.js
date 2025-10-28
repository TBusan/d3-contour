// Enhanced utility functions for d3-contour
export function slice(array, i) {
  return Array.prototype.slice.call(array, i);
}

export function ascending(a, b) {
  return a - b;
}

// Enhanced area calculation with better numerical stability
export function area(ring) {
  let i = 0, n = ring.length, sum = 0;
  let x0 = ring[n - 1][0];
  let y0 = ring[n - 1][1];
  let x1, y1;
  
  while (i < n) {
    x1 = ring[i][0];
    y1 = ring[i][1];
    sum += y0 * x1 - x0 * y1;
    x0 = x1;
    y0 = y1;
    i++;
  }
  
  return sum / 2;
}

export function constant(x) {
  return function() {
    return x;
  };
}

// Optimized spatial index for hole assignment
export class SpatialIndex {
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
    
    // First pass: bounding box test (fast rejection)
    for (let i = 0; i < this.bounds.length; i++) {
      const bounds = this.bounds[i];
      if (holeBounds.minX >= bounds.minX && holeBounds.maxX <= bounds.maxX &&
          holeBounds.minY >= bounds.minY && holeBounds.maxY <= bounds.maxY) {
        candidates.push(i);
      }
    }
    
    // Sort candidates by area (smallest containing polygon first)
    candidates.sort((a, b) => {
      const areaA = (this.bounds[a].maxX - this.bounds[a].minX) * 
                    (this.bounds[a].maxY - this.bounds[a].minY);
      const areaB = (this.bounds[b].maxX - this.bounds[b].minX) * 
                    (this.bounds[b].maxY - this.bounds[b].minY);
      return areaA - areaB;
    });
    
    // Second pass: actual containment test
    for (const idx of candidates) {
      if (contains(this.polygons[idx].polygon[0], hole) !== -1) {
        return idx;
      }
    }
    
    return -1;
  }
}

// Enhanced point-in-polygon test with optimizations
export function contains(ring, hole) {
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

export function noop() {}

export function finite(x) {
  return isFinite(x) ? x : NaN;
}

export function above(x, value) {
  return x != null && isFinite(x) && x >= value;
}

export function valid(v) {
  return v == null || isNaN(v = +v) ? -Infinity : v;
}

// Enhanced smoothing function with configurable intensity
export function createSmoothingFunction(intensity = 1.0) {
  if (intensity <= 0) return noop;
  
  return function smoothEnhanced(ring, values, value, dx, dy) {
    ring.forEach(function(point) {
      var x = point[0],
          y = point[1],
          xt = x | 0,
          yt = y | 0,
          v1 = valid(values[yt * dx + xt]);
      
      if (x > 0 && x < dx && xt === x) {
        point[0] = smooth1(x, valid(values[yt * dx + xt - 1]), v1, value, intensity);
      }
      if (y > 0 && y < dy && yt === y) {
        point[1] = smooth1(y, valid(values[(yt - 1) * dx + xt]), v1, value, intensity);
      }
    });
  };
}

function smooth1(x, v0, v1, value, intensity) {
  const a = value - v0;
  const b = v1 - v0;
  const d = isFinite(a) && isFinite(b) ? a / b : 0.5;
  
  // Apply intensity factor to the interpolation
  const t = Math.max(0, Math.min(1, 0.5 + (d - 0.5) * intensity));
  return x + t - 0.5;
}

// Saddle point disambiguation
export function disambiguateSaddle(caseIndex, corners, value) {
  if (caseIndex !== 5 && caseIndex !== 10) return caseIndex;
  
  // Calculate average value at center of the cell
  const avg = (corners[0] + corners[1] + corners[2] + corners[3]) / 4;
  
  // Determine which configuration to use based on the center value
  if (caseIndex === 5) {
    return value >= avg ? 5 : 10;
  } else { // caseIndex === 10
    return value >= avg ? 10 : 5;
  }
}

// Function to detect and fix self-intersections
export function hasSelfIntersections(ring) {
  const n = ring.length;
  if (n < 4) return false;
  
  // Check for self-intersections between non-adjacent segments
  for (let i = 0; i < n; i++) {
    const a = ring[i];
    const b = ring[(i + 1) % n];
    
    for (let j = i + 2; j < n; j++) {
      // Skip adjacent segments and last-first segment check if already checking first-second
      if (i === 0 && j === n - 1) continue;
      
      const c = ring[j];
      const d = ring[(j + 1) % n];
      
      if (segmentIntersection(a, b, c, d)) {
        return true;
      }
    }
  }
  return false;
}

// Check if two line segments intersect with more precise floating point handling
function segmentIntersection(a, b, c, d) {
  // Check if segments share an endpoint (already handled in existing implementation)
  if (pointsEqual(a, c) || pointsEqual(a, d) || pointsEqual(b, c) || pointsEqual(b, d)) {
    return false;
  }
  
  // Use robust orientation test for better numerical stability
  function orientation(p, q, r) {
    const val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1]);
    const epsilon = 1e-9; // Tolerance for floating point errors
    
    if (Math.abs(val) < epsilon) return 0;  // Collinear
    return (val > 0) ? 1 : 2; // Clockwise or counterclockwise
  }
  
  // Check if point q is on segment pr (assuming they're collinear)
  function onSegment(p, q, r) {
    return q[0] <= Math.max(p[0], r[0]) && q[0] >= Math.min(p[0], r[0]) &&
           q[1] <= Math.max(p[1], r[1]) && q[1] >= Math.min(p[1], r[1]);
  }
  
  // Find the orientation of triplets
  const o1 = orientation(a, b, c);
  const o2 = orientation(a, b, d);
  const o3 = orientation(c, d, a);
  const o4 = orientation(c, d, b);
  
  // General case of intersection
  if (o1 !== o2 && o3 !== o4) return true;
  
  // Special Cases: Collinear points
  if (o1 === 0 && onSegment(a, c, b)) return true;
  if (o2 === 0 && onSegment(a, d, b)) return true;
  if (o3 === 0 && onSegment(c, a, d)) return true;
  if (o4 === 0 && onSegment(c, b, d)) return true;
  
  return false;
}

function pointsEqual(p1, p2) {
  const epsilon = 1e-9;
  return Math.abs(p1[0] - p2[0]) < epsilon && Math.abs(p1[1] - p2[1]) < epsilon;
}

// Helper function to convert contours to GeoJSON
export function toGeoJSON(contourResult, properties = {}) {
  if (Array.isArray(contourResult)) {
    return {
      type: "FeatureCollection",
      features: contourResult.map((item, i) => toGeoJSON(item, { ...properties, index: i }))
    };
  }
  
  if (contourResult.type === "MultiLineString" || contourResult.type === "MultiPolygon") {
    return {
      type: "Feature",
      properties: { ...properties, value: contourResult.value },
      geometry: {
        type: contourResult.type,
        coordinates: contourResult.coordinates
      }
    };
  }
  
  // Handle 'both' mode result
  if (contourResult.lines && contourResult.surfaces) {
    return {
      type: "FeatureCollection",
      features: [
        {
          type: "Feature",
          properties: { ...properties, type: "lines", value: contourResult.value },
          geometry: {
            type: contourResult.lines.type,
            coordinates: contourResult.lines.coordinates
          }
        },
        {
          type: "Feature",
          properties: { ...properties, type: "surfaces", value: contourResult.value },
          geometry: {
            type: contourResult.surfaces.type,
            coordinates: contourResult.surfaces.coordinates
          }
        }
      ]
    };
  }
  
  return contourResult;
} 