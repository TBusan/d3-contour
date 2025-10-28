/**
 * Topology utility functions for contour generation
 * Handles self-intersections and other topological errors
 */

import { segmentIntersect } from './geometry.js';
import { pointsEqual } from './helpers.js';

// Check if a ring has self-intersections
export function hasSelfIntersections(ring) {
    if (ring.length < 4) return false;
    
    // Check each pair of non-adjacent line segments for intersection
    for (let i = 0; i < ring.length - 1; i++) {
        const a = ring[i];
        const b = ring[i + 1];
        
        for (let j = i + 2; j < ring.length - 1; j++) {
            // Skip adjacent segments
            if (j === i - 1 || j === i || j === i + 1) continue;
            
            const c = ring[j];
            const d = ring[j + 1];
            
            // Skip if the segments share an endpoint
            if (pointsEqual(a, c) || pointsEqual(a, d) || pointsEqual(b, c) || pointsEqual(b, d)) continue;
            
            // Check for intersection
            if (segmentIntersect(a, b, c, d)) {
                return true;
            }
        }
    }
    
    return false;
}

// Remove self-intersections from a line
export function removeSelfIntersections(line) {
    if (line.length < 4) return line;
    
    const result = [line[0]];
    let currentPoint = line[0];
    
    // Walk through the line, skipping points that would create self-intersections
    for (let i = 1; i < line.length; i++) {
        const nextPoint = line[i];
        let hasIntersection = false;
        
        // Check if adding this segment would create an intersection
        for (let j = 0; j < result.length - 1; j++) {
            if (segmentIntersect(
                currentPoint, nextPoint,
                result[j], result[j + 1]
            )) {
                hasIntersection = true;
                break;
            }
        }
        
        if (!hasIntersection) {
            result.push(nextPoint);
            currentPoint = nextPoint;
        }
    }
    
    return result;
}

// Simplify a polyline using the Douglas-Peucker algorithm
export function simplifyPolyline(points, tolerance = 1.0) {
    if (points.length <= 2) return points;
    
    // Find the point with the maximum distance
    let maxDistance = 0;
    let maxIndex = 0;
    
    const firstPoint = points[0];
    const lastPoint = points[points.length - 1];
    
    for (let i = 1; i < points.length - 1; i++) {
        const distance = perpendicularDistance(points[i], firstPoint, lastPoint);
        if (distance > maxDistance) {
            maxDistance = distance;
            maxIndex = i;
        }
    }
    
    // If max distance is greater than tolerance, recursively simplify
    if (maxDistance > tolerance) {
        const firstHalf = simplifyPolyline(points.slice(0, maxIndex + 1), tolerance);
        const secondHalf = simplifyPolyline(points.slice(maxIndex), tolerance);
        
        // Concatenate the two halves, avoiding duplicating the middle point
        return firstHalf.slice(0, -1).concat(secondHalf);
    }
    
    // Otherwise, return just the endpoints
    return [firstPoint, lastPoint];
}

// Calculate perpendicular distance from a point to a line segment
function perpendicularDistance(point, lineStart, lineEnd) {
    const [x, y] = point;
    const [x1, y1] = lineStart;
    const [x2, y2] = lineEnd;
    
    const dx = x2 - x1;
    const dy = y2 - y1;
    
    // If the line is just a point, return distance to that point
    const length = Math.sqrt(dx * dx + dy * dy);
    if (length === 0) return Math.sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));
    
    // Calculate the perpendicular distance
    const t = ((x - x1) * dx + (y - y1) * dy) / (length * length);
    
    if (t < 0) {
        // Beyond lineStart
        return Math.sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));
    }
    if (t > 1) {
        // Beyond lineEnd
        return Math.sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2));
    }
    
    // Perpendicular point on line
    const projX = x1 + t * dx;
    const projY = y1 + t * dy;
    
    return Math.sqrt((x - projX) * (x - projX) + (y - projY) * (y - projY));
} 