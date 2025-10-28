/**
 * Optimized spatial index for polygon and hole matching
 * Addresses the O(n²) performance issue in the original implementation
 */

import { contains } from './geometry.js';

export class SpatialIndex {
    constructor() {
        this.polygons = [];
        this.bounds = [];
    }
    
    // Add a polygon to the index with its bounding box
    addPolygon(polygon, index) {
        this.polygons.push({ polygon, index });
        this.bounds.push(this.computeBounds(polygon[0])); // Outer ring
    }
    
    // Calculate the bounding box of a ring
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
    
    // Find the smallest polygon that contains the hole
    findContainingPolygon(hole) {
        const holeBounds = this.computeBounds(hole);
        const candidates = [];
        
        // First pass: bounding box test (much faster than full containment test)
        for (let i = 0; i < this.bounds.length; i++) {
            const bounds = this.bounds[i];
            if (holeBounds.minX >= bounds.minX && holeBounds.maxX <= bounds.maxX &&
                holeBounds.minY >= bounds.minY && holeBounds.maxY <= bounds.maxY) {
                candidates.push(i);
            }
        }
        
        // If we have no candidates after the bounding box test, return early
        if (candidates.length === 0) {
            return -1;
        }
        
        // For multiple candidates, we need to find the smallest containing polygon
        // This handles nested polygons correctly by preferring the innermost container
        let bestIndex = -1;
        let smallestArea = Infinity;
        
        for (const idx of candidates) {
            const polygon = this.polygons[idx].polygon;
            const ring = polygon[0];
            
            // Check if the hole is contained within this polygon
            if (contains(ring, hole) !== -1) {
                // Calculate approximate area to determine the smallest container
                const bounds = this.bounds[idx];
                const area = (bounds.maxX - bounds.minX) * (bounds.maxY - bounds.minY);
                
                if (area < smallestArea) {
                    smallestArea = area;
                    bestIndex = this.polygons[idx].index;
                }
            }
        }
        
        return bestIndex;
    }
} 