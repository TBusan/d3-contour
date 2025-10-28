// Data adapter to handle different input formats including testData1.js and testData2.js
export default function() {
  
  function adaptData(data) {
    // Handle structured data format like testData1.js and testData2.js
    if (data && typeof data === 'object' && data.x && data.y && data.v) {
      return adaptStructuredData(data);
    }
    
    // Handle flat array format
    if (Array.isArray(data)) {
      return {
        values: data,
        width: Math.sqrt(data.length),
        height: Math.sqrt(data.length)
      };
    }
    
    throw new Error("Unsupported data format");
  }
  
  function adaptStructuredData(data) {
    const {x, y, v} = data;
    const width = x.length;
    const height = y.length;
    
    // Flatten the 2D v array to 1D array in row-major order
    const values = new Array(width * height);
    
    for (let j = 0; j < height; j++) {
      for (let i = 0; i < width; i++) {
        values[j * width + i] = v[j][i];
      }
    }
    
    return {
      values: values,
      width: width,
      height: height,
      x: x,
      y: y,
      xmin: x[0],
      xmax: x[x.length - 1],
      ymin: y[0],
      ymax: y[y.length - 1],
      dx: width > 1 ? (x[x.length - 1] - x[0]) / (width - 1) : 1,
      dy: height > 1 ? (y[y.length - 1] - y[0]) / (height - 1) : 1
    };
  }
  
  function createContourData(adaptedData, contourResult) {
    // Transform contour coordinates back to original coordinate system
    if (adaptedData.x && adaptedData.y) {
      return transformCoordinates(contourResult, adaptedData);
    }
    return contourResult;
  }
  
  function transformCoordinates(contourResult, adaptedData) {
    const {x, y, width, height} = adaptedData;
    
    return contourResult.map(contour => ({
      ...contour,
      coordinates: contour.coordinates.map(polygon => 
        polygon.map(ring => 
          ring.map(point => [
            interpolateCoordinate(point[0], x, width),
            interpolateCoordinate(point[1], y, height)
          ])
        )
      )
    }));
  }
  
  function interpolateCoordinate(gridIndex, coords, size) {
    // Linear interpolation between grid coordinates
    const i = Math.floor(gridIndex);
    const f = gridIndex - i;
    
    if (i < 0) return coords[0];
    if (i >= size - 1) return coords[size - 1];
    
    return coords[i] * (1 - f) + coords[i + 1] * f;
  }
  
  adaptData.adaptStructuredData = adaptStructuredData;
  adaptData.createContourData = createContourData;
  adaptData.transformCoordinates = transformCoordinates;
  adaptData.interpolateCoordinate = interpolateCoordinate;
  
  return adaptData;
}