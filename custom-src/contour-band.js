import contours from "./contours.js";
import {extent, nice, thresholdSturges, ticks} from "d3-array";
import {slice} from "./array.js";
import ascending from "./ascending.js";
import constant from "./constant.js";

export default function() {
  var dx = 1,
      dy = 1,
      threshold = thresholdSturges,
      contourInstance = contours();
  
  function contourBand(values) {
    // Filter out null values for threshold calculation
    const validValues = values.filter(v => v != null && !isNaN(+v));
    if (validValues.length === 0) return []; // No valid data points
    
    var tz = threshold(validValues);

    // Convert number of thresholds into uniform thresholds
    if (!Array.isArray(tz)) {
      const e = extent(validValues, finite);
      if (e[0] === undefined || e[1] === undefined) {
        return []; // No valid data
      }
      tz = ticks(...nice(e[0], e[1], tz), tz);
      while (tz[tz.length - 1] >= e[1]) tz.pop();
      while (tz.length > 1 && tz[0] < e[0]) tz.shift();
    } else {
      tz = tz.slice().sort(ascending);
    }

    // Generate bands between adjacent threshold values
    const bands = [];
    for (let i = 0; i < tz.length - 1; i++) {
      const lowerValue = tz[i];
      const upperValue = tz[i + 1];
      
      const band = {
        type: "MultiPolygon",
        lowerValue: lowerValue,
        upperValue: upperValue,
        coordinates: generateBandGeometry(values, lowerValue, upperValue)
      };
      
      bands.push(band);
    }

    return bands;
  }

  function generateBandGeometry(values, lowerValue, upperValue) {
    // Get the contour polygons for both thresholds
    const lowerContour = contourInstance.contour(values, lowerValue);
    const upperContour = contourInstance.contour(values, upperValue);
    
    // For band generation, we need to invert the "hole" status of the upper contour
    // The upper contour becomes the inner boundary of the band
    const lowerPolygons = lowerContour.coordinates; 
    const upperPolygons = invertHoles(upperContour.coordinates);

    // Combine polygons to form bands
    return combineContours(lowerPolygons, upperPolygons);
  }

  // Invert the orientation of polygons - exterior becomes interior and vice versa
  function invertHoles(polygons) {
    return polygons.map(polygon => {
      // Reverse the order of rings (first becomes hole, holes become exteriors)
      if (polygon.length > 1) {
        return [polygon[0]].concat(polygon.slice(1).reverse());
      }
      return polygon;
    });
  }

  // Combine lower and upper contours to form bands
  function combineContours(lowerPolygons, upperPolygons) {
    const result = [];
    
    // In the simple case, we can use the lower contour exterior with the upper contour holes
    for (const lowerPoly of lowerPolygons) {
      const exterior = lowerPoly[0]; // Exterior ring of lower polygon
      const holes = [];
      
      // Find upper polygon rings that are inside this lower polygon exterior
      for (const upperPoly of upperPolygons) {
        for (const ring of upperPoly) {
          // Check if the upper ring is inside the lower exterior
          if (isRingInside(ring, exterior)) {
            holes.push(ring);
          }
        }
      }
      
      // Create a new polygon with the lower exterior and upper holes
      result.push([exterior, ...holes]);
    }
    
    return result;
  }

  // Determine if ring1 is inside ring2 using a simple point-in-polygon test
  function isRingInside(ring1, ring2) {
    // Use the first point of ring1 to test
    const point = ring1[0];
    return isPointInPolygon(point, ring2);
  }

  // Point-in-polygon test using ray casting algorithm
  function isPointInPolygon(point, polygon) {
    const x = point[0], y = point[1];
    let inside = false;
    
    for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
      const xi = polygon[i][0], yi = polygon[i][1];
      const xj = polygon[j][0], yj = polygon[j][1];
      
      const intersect = ((yi > y) !== (yj > y))
          && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
      if (intersect) inside = !inside;
    }
    
    return inside;
  }

  // When computing the extent, ignore invalid values
  function finite(x) {
    return isFinite(x) ? x : NaN;
  }

  contourBand.size = function(_) {
    if (!arguments.length) return [dx, dy];
    var _0 = Math.floor(_[0]), _1 = Math.floor(_[1]);
    if (!(_0 >= 0 && _1 >= 0)) throw new Error("invalid size");
    contourInstance.size(_);
    return dx = _0, dy = _1, contourBand;
  };

  contourBand.thresholds = function(_) {
    return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? constant(slice.call(_)) : constant(_), contourBand) : threshold;
  };

  contourBand.smooth = function(_) {
    if (!arguments.length) return contourInstance.smooth();
    contourInstance.smooth(_);
    return contourBand;
  };

  return contourBand;
} 