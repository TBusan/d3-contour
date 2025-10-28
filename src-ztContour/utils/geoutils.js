/**
 * GeoJSON export utilities for contour generation
 * Provides tools to convert contour data to standard GeoJSON formats
 */

// Convert contour output to GeoJSON feature or feature collection
// Works with lines, surfaces or both
export default function toGeoJSON(contourResult, properties = {}) {
    if (Array.isArray(contourResult)) {
        return {
            type: "FeatureCollection",
            features: contourResult.map((item, i) => toGeoJSON(item, { ...properties, index: i }))
        };
    }
    
    if (contourResult.lines && contourResult.surfaces) {
        // Both mode
        return {
            type: "FeatureCollection",
            features: [
                {
                    type: "Feature",
                    properties: { ...properties, type: "lines", value: contourResult.value },
                    geometry: contourResult.lines
                },
                {
                    type: "Feature",
                    properties: { ...properties, type: "surfaces", value: contourResult.value },
                    geometry: contourResult.surfaces
                }
            ]
        };
    }
    
    return {
        type: "Feature",
        properties: { ...properties, value: contourResult.value },
        geometry: contourResult
    };
}

// Apply a projection function to contour geometry
export function transformContours(contourResult, projectionFn) {
    if (Array.isArray(contourResult)) {
        return contourResult.map(item => transformContours(item, projectionFn));
    }
    
    if (contourResult.lines && contourResult.surfaces) {
        return {
            lines: transformGeometry(contourResult.lines, projectionFn),
            surfaces: transformGeometry(contourResult.surfaces, projectionFn),
            value: contourResult.value
        };
    }
    
    return transformGeometry(contourResult, projectionFn);
}

// Transform a GeoJSON geometry with a projection function
export function transformGeometry(geometry, projectionFn) {
    if (!geometry || !geometry.coordinates) return geometry;
    
    const transformed = {
        ...geometry,
        coordinates: transformCoordinates(geometry.coordinates, projectionFn)
    };
    
    return transformed;
}

// Transform coordinates array recursively
function transformCoordinates(coordinates, projectionFn) {
    if (!Array.isArray(coordinates)) return coordinates;
    
    if (Array.isArray(coordinates[0])) {
        if (typeof coordinates[0][0] === 'number' && coordinates[0].length >= 2) {
            // This is a point, apply projection
            return projectionFn(coordinates);
        }
        // Recurse into nested arrays
        return coordinates.map(coord => transformCoordinates(coord, projectionFn));
    }
    
    return coordinates;
} 