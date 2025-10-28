import {constant} from "./utils.js";
import contours from "./contours.js";

export default function() {
  var x = d => d[0],
      y = d => d[1],
      weight = constant(1),
      dx = 960,
      dy = 500,
      thresholds = 20,
      bandwidth = 20.4939015319192,
      kernelSize = 3,
      kernelDensityFactor = 1.0,
      nullValuesMask = false,
      mode = "surfaces";

  function density(data) {
    var values = new Float64Array(dx * dy),
        nullMask = nullValuesMask ? new Uint8Array(dx * dy) : null,
        i = -1,
        n = data.length,
        xi,
        yi,
        wi;

    // Initialize grid to zeros
    while (++i < dx * dy) {
      values[i] = 0;
      if (nullMask) nullMask[i] = 1; // Initially all values are valid
    }

    // Compute kernel density estimation (KDE)
    i = -1;
    while (++i < n) {
      xi = x(data[i]);
      yi = y(data[i]);
      wi = weight(data[i]);
      
      // Skip invalid data points
      if (isNaN(xi) || isNaN(yi) || isNaN(wi)) continue;
      
      // Convert to grid coordinates
      let gx = (xi - bandwidth / 2) / bandwidth;
      let gy = (yi - bandwidth / 2) / bandwidth;
      
      // Floor to get grid cell
      let x0 = Math.floor(gx);
      let y0 = Math.floor(gy);
      
      // Apply kernel density estimation with given kernel size
      applyKernel(values, x0, y0, gx - x0, gy - y0, wi * kernelDensityFactor, dx, dy);
    }

    // Apply threshold function and generate contours
    const contourGenerator = contours()
        .size([dx, dy])
        .thresholds(thresholds)
        .mode(mode);
    
    // Add null mask if enabled
    if (nullValuesMask) {
      // Mark cells with zero density as null values
      for (let i = 0; i < values.length; i++) {
        if (values[i] <= 0) nullMask[i] = 0;
      }
      contourGenerator.nullMask(nullMask);
    }
    
    return contourGenerator(values);
  }
  
  // Apply a more accurate kernel for density estimation
  function applyKernel(values, x0, y0, xf, yf, weight, gridDx, gridDy) {
    const kernelRadius = Math.floor(kernelSize / 2);
    const sigma = kernelSize / 6; // Ensure 3-sigma rule covers the kernel
    
    // Calculate normalization factor
    let normalization = 0;
    
    // For each point in the kernel
    for (let ky = -kernelRadius; ky <= kernelRadius; ky++) {
      for (let kx = -kernelRadius; kx <= kernelRadius; kx++) {
        // Calculate grid coordinates
        const gx = x0 + kx;
        const gy = y0 + ky;
        
        // Skip out-of-bounds cells
        if (gx < 0 || gx >= gridDx || gy < 0 || gy >= gridDy) continue;
        
        // Calculate normalized distance from kernel center
        const distX = (kx - xf);
        const distY = (ky - yf);
        const distance = Math.sqrt(distX * distX + distY * distY);
        
        // Apply Gaussian kernel
        const kernelValue = Math.exp(-(distance * distance) / (2 * sigma * sigma));
        normalization += kernelValue;
        
        // Apply to grid
        values[gy * gridDx + gx] += kernelValue * weight;
      }
    }
    
    // Normalize if needed
    if (normalization > 0 && kernelDensityFactor > 0) {
      // Apply post-normalization to ensure consistent weighting
      const normalizedWeight = weight / normalization;
      
      for (let ky = -kernelRadius; ky <= kernelRadius; ky++) {
        for (let kx = -kernelRadius; kx <= kernelRadius; kx++) {
          const gx = x0 + kx;
          const gy = y0 + ky;
          
          if (gx < 0 || gx >= gridDx || gy < 0 || gy >= gridDy) continue;
          
          const distX = (kx - xf);
          const distY = (ky - yf);
          const distance = Math.sqrt(distX * distX + distY * distY);
          
          const kernelValue = Math.exp(-(distance * distance) / (2 * sigma * sigma));
          values[gy * gridDx + gx] += (kernelValue * normalizedWeight) - (kernelValue * weight);
        }
      }
    }
  }

  density.x = function(_) {
    return arguments.length ? (x = typeof _ === "function" ? _ : constant(+_), density) : x;
  };

  density.y = function(_) {
    return arguments.length ? (y = typeof _ === "function" ? _ : constant(+_), density) : y;
  };

  density.weight = function(_) {
    return arguments.length ? (weight = typeof _ === "function" ? _ : constant(+_), density) : weight;
  };

  density.size = function(_) {
    if (!arguments.length) return [dx, dy];
    var _0 = +_[0], _1 = +_[1];
    if (!(_0 >= 0 && _1 >= 0)) throw new Error("invalid size");
    return dx = _0, dy = _1, density;
  };

  density.cellSize = function(_) {
    if (!arguments.length) return bandwidth;
    if (!((_ = +_) >= 0)) throw new Error("invalid cell size");
    return bandwidth = _, density;
  };

  density.thresholds = function(_) {
    return arguments.length ? (thresholds = typeof _ === "function" ? _ : Array.isArray(_) ? constant(_.slice()) : constant(+_), density) : thresholds;
  };

  density.bandwidth = function(_) {
    if (!arguments.length) return Math.sqrt(bandwidth * bandwidth * 20);
    if (!((_ = +_) >= 0)) throw new Error("invalid bandwidth");
    return bandwidth = _ / Math.sqrt(20), density;
  };
  
  density.kernelSize = function(_) {
    if (!arguments.length) return kernelSize;
    if (!((_ = +_) >= 1)) throw new Error("invalid kernel size");
    return kernelSize = Math.floor(_), density;
  };
  
  density.kernelDensityFactor = function(_) {
    if (!arguments.length) return kernelDensityFactor;
    if (!((_ = +_) >= 0)) throw new Error("invalid kernel density factor");
    return kernelDensityFactor = _, density;
  };
  
  density.nullValuesMask = function(_) {
    if (!arguments.length) return nullValuesMask;
    nullValuesMask = !!_;
    return density;
  };
  
  density.mode = function(_) {
    if (!arguments.length) return mode;
    mode = _ + "";
    return density;
  };

  return density;
} 