/**
 * Helper functions for the enhanced contour library
 */

// No-operation function
export function noop() {}

// When computing the extent, ignore infinite values (as well as invalid ones).
export function finite(x) {
    return isFinite(x) ? x : NaN;
}

// Is the (possibly invalid) x greater than or equal to the (known valid) value?
// Treat any invalid value as below negative infinity.
export function above(x, value) {
    return x != null && isFinite(x) && x >= value;
}

// During smoothing, treat any invalid value as negative infinity.
export function valid(v) {
    return v == null || isNaN(v = +v) ? -Infinity : v;
}

// Calculate distance between two points
export function distance(p1, p2) {
    const dx = p1[0] - p2[0];
    const dy = p1[1] - p2[1];
    return Math.sqrt(dx * dx + dy * dy);
}

// Check if two points are equal (within a small epsilon)
export function pointsEqual(p1, p2) {
    const epsilon = 1e-6;
    return Math.abs(p1[0] - p2[0]) < epsilon && Math.abs(p1[1] - p2[1]) < epsilon;
} 