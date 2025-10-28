/**
 * Geometric utility functions for contour generation
 */

// Calculate the area of a ring (polygon)
// Positive for counter-clockwise (exterior rings), negative for clockwise (hole rings)
export function area(ring) {
    let i = 0, n = ring.length;
    let area = ring[n - 1][1] * ring[0][0] - ring[n - 1][0] * ring[0][1];
    while (++i < n) area += ring[i - 1][1] * ring[i][0] - ring[i - 1][0] * ring[i][1];
    return area;
}

// Check if a ring contains a point
// Returns -1 for contained points, 0 for boundary points, 1 for exterior points
export function contains(ring, hole) {
    let i = -1, n = hole.length, c;
    while (++i < n) if (c = ringContains(ring, hole[i])) return c;
    return 0;
}

// Helper function for contains() - tests if a single point is inside a ring
export function ringContains(ring, point) {
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

// Check if a point is on a line segment
export function segmentContains(a, b, c) {
    const i = +(a[0] === b[0]);
    return collinear(a, b, c) && within(a[i], c[i], b[i]);
}

// Check if three points are collinear (on the same line)
export function collinear(a, b, c) {
    return (b[0] - a[0]) * (c[1] - a[1]) === (c[0] - a[0]) * (b[1] - a[1]);
}

// Check if a value is between two other values
export function within(p, q, r) {
    return p <= q && q <= r || r <= q && q <= p;
}

// Check if two line segments intersect
export function segmentIntersect(a, b, c, d) {
    // Calculate the direction vectors
    const ab = [b[0] - a[0], b[1] - a[1]];
    const cd = [d[0] - c[0], d[1] - c[1]];
    
    // Calculate the determinant
    const det = ab[0] * cd[1] - ab[1] * cd[0];
    
    // If determinant is zero, lines are parallel
    if (Math.abs(det) < 1e-10) return false;
    
    // Calculate the parameters t and s
    const ac = [c[0] - a[0], c[1] - a[1]];
    const t = (ac[0] * cd[1] - ac[1] * cd[0]) / det;
    const s = (ac[0] * ab[1] - ac[1] * ab[0]) / det;
    
    // Intersection occurs if t and s are both in [0,1]
    return t >= 0 && t <= 1 && s >= 0 && s <= 1;
} 