/**
 * Enhanced D3-Contour Library
 * Improved version addressing the limitations of the original d3-contour
 * 
 * Features:
 * - Support for both contour lines and contour surfaces
 * - Saddle point disambiguation
 * - Optimized hole assignment algorithm
 * - Configurable smoothing
 * - Null value handling and boundary clipping
 * - GeoJSON export for both lines and polygons
 */

export { default as contours } from "./contours.js";
export { default as contourDensity } from "./density.js";
export { default as toGeoJSON } from "./utils/geoutils.js"; 