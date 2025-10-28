/**
 * Enhanced smoothing functions for contour generation
 * Provides configurable smoothing with multiple strategies
 */

import { valid } from './helpers.js';

// Create a smoothing function with configurable parameters
export function createSmoothing(smoothingLevel = 1.0) {
    if (smoothingLevel <= 0) return function() {}; // No smoothing
    
    return function smoothAdvanced(ring, values, value, dx, dy) {
        // Threshold for detecting grid-aligned points that need smoothing
        const SMOOTH_THRESHOLD = 0.03 * smoothingLevel; 
        
        ring.forEach(function(point) {
            const x = point[0];
            const y = point[1];
            const xt = Math.floor(x);
            const yt = Math.floor(y);
            
            if (xt >= 0 && xt < dx - 1 && yt >= 0 && yt < dy - 1) {
                // Get fractional position within grid cell
                const xFrac = x - xt;
                const yFrac = y - yt;
                
                // Detect and smooth grid-aligned points to avoid artifacts
                if (Math.abs(xFrac) < SMOOTH_THRESHOLD || Math.abs(xFrac - 1) < SMOOTH_THRESHOLD) {
                    // Smoothly adjust x coordinate away from grid lines
                    point[0] += (xFrac < 0.5 ? 1 : -1) * SMOOTH_THRESHOLD * smoothingLevel * 0.5;
                }
                
                if (Math.abs(yFrac) < SMOOTH_THRESHOLD || Math.abs(yFrac - 1) < SMOOTH_THRESHOLD) {
                    // Smoothly adjust y coordinate away from grid lines
                    point[1] += (yFrac < 0.5 ? 1 : -1) * SMOOTH_THRESHOLD * smoothingLevel * 0.5;
                }
                
                // Enhanced interpolation at cell edges
                if (Math.abs(xFrac - 0.5) < SMOOTH_THRESHOLD && x > 0 && x < dx) {
                    const v0 = valid(values[yt * dx + xt]);
                    const v1 = valid(values[yt * dx + xt + 1]);
                    point[0] = smoothInterpolate(x, v0, v1, value, smoothingLevel);
                }
                
                if (Math.abs(yFrac - 0.5) < SMOOTH_THRESHOLD && y > 0 && y < dy) {
                    const v0 = valid(values[yt * dx + xt]);
                    const v1 = valid(values[(yt + 1) * dx + xt]);
                    point[1] = smoothInterpolate(y, v0, v1, value, smoothingLevel);
                }
            }
        });
    };
}

// Enhanced interpolation function with smoothing factor
function smoothInterpolate(coord, v0, v1, value, smoothingFactor = 1.0) {
    // Get the base coordinate (integer part)
    const base = Math.floor(coord);
    
    // Handle edge cases
    if (v0 === v1) return base + 0.5;
    if (!isFinite(v0) || !isFinite(v1)) return coord;
    
    // Calculate interpolation factor
    const a = value - v0;
    const b = v1 - v0;
    let d = isFinite(a) && isFinite(b) && b !== 0 ? a / b : 0.5;
    
    // Apply smoothing factor - higher values make interpolation more aggressive
    if (smoothingFactor !== 1.0) {
        // Move interpolation factor toward 0.5 for smoother transitions
        d = d + (0.5 - d) * (1 - smoothingFactor);
    }
    
    // Keep within valid range
    if (d < 0) d = 0;
    if (d > 1) d = 1;
    
    // Apply simple linear interpolation
    return base + d;
} 