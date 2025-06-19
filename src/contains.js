export default function(ring, point) {
  // Add safety check for inputs
  if (!Array.isArray(ring) || !Array.isArray(point) || point.length < 2) {
    return -1;
  }

  const n = ring.length;
  let p = ring[0];
  let x = point[0], y = point[1];
  let x0 = p[0], y0 = p[1];
  let x1, y1;
  let inside = -1;
  for (let i = 1; i <= n; ++i) {
    p = ring[i % n];
    // Add safety check for each point in the ring
    if (!Array.isArray(p) || p.length < 2) continue;
    
    x1 = p[0];
    y1 = p[1];
    if (((y1 > y) !== (y0 > y)) && ((x0 + (y - y0) * (x1 - x0) / (y1 - y0)) > x)) inside = -inside;
    x0 = x1, y0 = y1;
  }
  return inside;
}

// Calculate the bounds of a polygon
function getBounds(points) {
  let minX = Infinity, minY = Infinity;
  let maxX = -Infinity, maxY = -Infinity;
  
  for (let i = 0; i < points.length; i++) {
    const [x, y] = points[i];
    if (x < minX) minX = x;
    if (y < minY) minY = y;
    if (x > maxX) maxX = x;
    if (y > maxY) maxY = y;
  }
  
  return [minX, minY, maxX, maxY];
}

// Check if two bounding boxes overlap
function boundsOverlap(bounds1, bounds2) {
  const [minX1, minY1, maxX1, maxY1] = bounds1;
  const [minX2, minY2, maxX2, maxY2] = bounds2;
  
  return !(
    maxX1 < minX2 ||
    maxX2 < minX1 ||
    maxY1 < minY2 ||
    maxY2 < minY1
  );
}

// Optimized point-in-polygon test using ray casting algorithm
function ringContains(ring, point) {
  var x = point[0], y = point[1], contains = -1;
  
  // Quick check if point is on any segment
  for (var i = 0, n = ring.length, j = n - 1; i < n; j = i++) {
    var pi = ring[i], xi = pi[0], yi = pi[1], 
        pj = ring[j], xj = pj[0], yj = pj[1];
    
    // If point is on the boundary, return 0
    if (segmentContains(pi, pj, point)) return 0;
    
    // Ray casting algorithm: count crossings
    // Only count crossings where the y-coordinate of the point is between segment y-coordinates
    // and the x-coordinate of the point is to the left of the intersection
    if (((yi > y) !== (yj > y)) && ((x < (xj - xi) * (y - yi) / (yj - yi) + xi))) {
      contains = -contains;
    }
  }
  
  return contains;
}

// Check if point lies on a line segment
function segmentContains(a, b, c) {
  var i = +(a[0] === b[0]); // Choose coordinate to check based on segment orientation
  return collinear(a, b, c) && within(a[i], c[i], b[i]);
}

// Check if three points are collinear (lie on the same line)
function collinear(a, b, c) {
  // Cross product will be 0 if points are collinear
  return Math.abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) < 1e-10;
}

// Check if a value is between two other values (inclusive)
function within(p, q, r) {
  return p <= q && q <= r || r <= q && q <= p;
}
