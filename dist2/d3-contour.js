// https://d3js.org/d3-contour/ v4.0.2 Copyright 2012-2023 Mike Bostock
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports, require('d3-array')) :
typeof define === 'function' && define.amd ? define(['exports', 'd3-array'], factory) :
(global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.d3 = global.d3 || {}, global.d3));
})(this, (function (exports, d3Array) { 'use strict';

var array = Array.prototype;

var slice = array.slice;

function ascending(a, b) {
  return a - b;
}

function area(ring) {
  var i = 0, n = ring.length, area = ring[n - 1][1] * ring[0][0] - ring[n - 1][0] * ring[0][1];
  while (++i < n) area += ring[i - 1][1] * ring[i][0] - ring[i - 1][0] * ring[i][1];
  return area;
}

var constant = x => () => x;

function contains(ring, hole) {
  var i = -1, n = hole.length, c;
  while (++i < n) if (c = ringContains(ring, hole[i])) return c;
  return 0;
}

function ringContains(ring, point) {
  var x = point[0], y = point[1], contains = -1;
  for (var i = 0, n = ring.length, j = n - 1; i < n; j = i++) {
    var pi = ring[i], xi = pi[0], yi = pi[1], pj = ring[j], xj = pj[0], yj = pj[1];
    if (segmentContains(pi, pj, point)) return 0;
    if (((yi > y) !== (yj > y)) && ((x < (xj - xi) * (y - yi) / (yj - yi) + xi))) contains = -contains;
  }
  return contains;
}

function segmentContains(a, b, c) {
  var i; return collinear(a, b, c) && within(a[i = +(a[0] === b[0])], c[i], b[i]);
}

function collinear(a, b, c) {
  return (b[0] - a[0]) * (c[1] - a[1]) === (c[0] - a[0]) * (b[1] - a[1]);
}

function within(p, q, r) {
  return p <= q && q <= r || r <= q && q <= p;
}

function noop() {}

// Marching squares case lookup table
// These define how to connect points to form contour lines for each possible corner configuration
var cases = [
  [], // Case 0: All corners are below threshold
  [[[1.0, 1.5], [0.5, 1.0]]], // Case 1
  [[[1.5, 1.0], [1.0, 1.5]]], // Case 2
  [[[1.5, 1.0], [0.5, 1.0]]], // Case 3
  [[[1.0, 0.5], [1.5, 1.0]]], // Case 4
  [[[1.0, 1.5], [0.5, 1.0]], [[1.0, 0.5], [1.5, 1.0]]], // Case 5: Saddle point
  [[[1.0, 0.5], [1.0, 1.5]]], // Case 6
  [[[1.0, 0.5], [0.5, 1.0]]], // Case 7
  [[[0.5, 1.0], [1.0, 0.5]]], // Case 8
  [[[1.0, 1.5], [1.0, 0.5]]], // Case 9
  [[[0.5, 1.0], [1.0, 0.5]], [[1.5, 1.0], [1.0, 1.5]]], // Case 10: Saddle point
  [[[1.5, 1.0], [1.0, 0.5]]], // Case 11
  [[[0.5, 1.0], [1.5, 1.0]]], // Case 12
  [[[1.0, 1.5], [1.5, 1.0]]], // Case 13
  [[[0.5, 1.0], [1.0, 1.5]]], // Case 14
  [] // Case 15: All corners are above threshold
];

function contours() {
  var dx = 1, // Width of each grid cell
      dy = 1, // Height of each grid cell
      threshold = d3Array.thresholdSturges,
      smooth = smoothLinear;

  function contours(values) {
    var tz = threshold(values);

    // Convert number of thresholds into uniform thresholds.
    if (!Array.isArray(tz)) {
      const validValues = values.filter(v => v != null && !isNaN(v));
      const e = d3Array.extent(validValues, finite);
      if (e[0] === undefined || e[1] === undefined) {
        return []; // No valid data points
      }
      tz = d3Array.ticks(...d3Array.nice(e[0], e[1], tz), tz);
      while (tz[tz.length - 1] >= e[1]) tz.pop();
      while (tz[1] < e[0]) tz.shift();
    } else {
      tz = tz.slice().sort(ascending);
    }

    return tz.map(value => contour(values, value));
  }

  // Accumulate, smooth contour rings, assign holes to exterior rings.
  function contour(values, value) {
    const v = value == null ? NaN : +value;
    if (isNaN(v)) throw new Error(`invalid value: ${value}`);

    var polygons = [],
        holes = [];

    isorings(values, v, function(ring) {
      if (ring.length > 2) { // Ignore degenerate rings
        smooth(ring, values, v);
        if (area(ring) > 0) polygons.push([ring]);
        else holes.push(ring);
      }
    });

    holes.forEach(function(hole) {
      for (var i = 0, n = polygons.length, polygon; i < n; ++i) {
        if (contains((polygon = polygons[i])[0], hole) !== -1) {
          polygon.push(hole);
          return;
        }
      }
    });

    return {
      type: "MultiPolygon",
      value: value,
      coordinates: polygons
    };
  }

  // Marching squares with isolines stitched into rings.
  function isorings(values, value, callback) {
    var fragmentByStart = new Array,
        fragmentByEnd = new Array,
        x, y, t0, t1, t2, t3;

    // Special case for the first row (y = -1, t2 = t3 = 0).
    x = y = -1;
    t1 = above(values[0], value);
    cases[t1 << 1].forEach(stitch);
    while (++x < dx - 1) {
      t0 = t1, t1 = above(values[x + 1], value);
      cases[t0 | t1 << 1].forEach(stitch);
    }
    cases[t1 << 0].forEach(stitch);

    // General case for the intermediate rows.
    while (++y < dy - 1) {
      x = -1;
      t1 = above(values[(y + 1) * dx], value);
      t2 = above(values[y * dx], value);
      cases[t1 << 1 | t2 << 2].forEach(stitch);
      while (++x < dx - 1) {
        t0 = t1, t1 = above(values[(y + 1) * dx + x + 1], value);
        t3 = t2, t2 = above(values[y * dx + x + 1], value);
        
        // Special handling for ambiguous cases (saddle points)
        const config = t0 | t1 << 1 | t2 << 2 | t3 << 3;
        if (config === 5 || config === 10) {
          // For saddle points, use the average of cell values to determine configuration
          const avg = averageValidValues([
            values[y * dx + x],
            values[y * dx + x + 1], 
            values[(y + 1) * dx + x], 
            values[(y + 1) * dx + x + 1]
          ]);
          
          // Use the average value compared to threshold to resolve ambiguity
          if ((config === 5 && avg < value) || (config === 10 && avg < value)) {
            // Original configuration
            cases[config].forEach(stitch);
          } else {
            // Alternate configuration
            cases[15 - config].forEach(stitch);
          }
        } else {
          // Normal cases
          cases[config].forEach(stitch);
        }
      }
      cases[t1 | t2 << 3].forEach(stitch);
    }

    // Special case for the last row (y = dy - 1, t0 = t1 = 0).
    x = -1;
    t2 = above(values[y * dx], value);
    cases[t2 << 2].forEach(stitch);
    while (++x < dx - 1) {
      t3 = t2, t2 = above(values[y * dx + x + 1], value);
      cases[t2 << 2 | t3 << 3].forEach(stitch);
    }
    cases[t2 << 3].forEach(stitch);

    function stitch(line) {
      var start = [line[0][0] + x, line[0][1] + y],
          end = [line[1][0] + x, line[1][1] + y],
          startIndex = index(start),
          endIndex = index(end),
          f, g;
      if (f = fragmentByEnd[startIndex]) {
        if (g = fragmentByStart[endIndex]) {
          delete fragmentByEnd[f.end];
          delete fragmentByStart[g.start];
          if (f === g) {
            f.ring.push(end);
            callback(f.ring);
          } else {
            fragmentByStart[f.start] = fragmentByEnd[g.end] = {start: f.start, end: g.end, ring: f.ring.concat(g.ring)};
          }
        } else {
          delete fragmentByEnd[f.end];
          f.ring.push(end);
          fragmentByEnd[f.end = endIndex] = f;
        }
      } else if (f = fragmentByStart[endIndex]) {
        if (g = fragmentByEnd[startIndex]) {
          delete fragmentByStart[f.start];
          delete fragmentByEnd[g.end];
          if (f === g) {
            f.ring.push(end);
            callback(f.ring);
          } else {
            fragmentByStart[g.start] = fragmentByEnd[f.end] = {start: g.start, end: f.end, ring: g.ring.concat(f.ring)};
          }
        } else {
          delete fragmentByStart[f.start];
          f.ring.unshift(start);
          fragmentByStart[f.start = startIndex] = f;
        }
      } else {
        fragmentByStart[startIndex] = fragmentByEnd[endIndex] = {start: startIndex, end: endIndex, ring: [start, end]};
      }
    }
  }

  function index(point) {
    return point[0] * 2 + point[1] * (dx + 1) * 4;
  }

  function averageValidValues(values) {
    const validValues = values.filter(v => v != null && !isNaN(+v));
    if (validValues.length === 0) return -Infinity; // Default below threshold if no valid values
    return validValues.reduce((sum, v) => sum + (+v), 0) / validValues.length;
  }

  // Enhanced smoothing that handles null values
  function smoothLinear(ring, values, value) {
    ring.forEach(function(point) {
      var x = point[0],
          y = point[1],
          xt = x | 0,
          yt = y | 0;
          
      // Ensure coordinates are within bounds
      if (xt < 0 || yt < 0 || xt >= dx || yt >= dy) return;
      
      const idx = yt * dx + xt;
      const v1 = valid(values[idx]);
      
      if (x > 0 && x < dx && xt === x) {
        const leftIdx = yt * dx + (xt - 1);
        const leftVal = valid(values[leftIdx]);
        
        if (isFinite(leftVal) && isFinite(v1)) {
          point[0] = smooth1(x, leftVal, v1, value);
        }
      }
      
      if (y > 0 && y < dy && yt === y) {
        const topIdx = (yt - 1) * dx + xt;
        const topVal = valid(values[topIdx]);
        
        if (isFinite(topVal) && isFinite(v1)) {
          point[1] = smooth1(y, topVal, v1, value);
        }
      }
    });
  }

  contours.contour = contour;

  contours.size = function(_) {
    if (!arguments.length) return [dx, dy];
    var _0 = Math.floor(_[0]), _1 = Math.floor(_[1]);
    if (!(_0 >= 0 && _1 >= 0)) throw new Error("invalid size");
    return dx = _0, dy = _1, contours;
  };

  contours.thresholds = function(_) {
    return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? constant(slice.call(_)) : constant(_), contours) : threshold;
  };

  contours.smooth = function(_) {
    return arguments.length ? (smooth = _ ? smoothLinear : noop, contours) : smooth === smoothLinear;
  };

  return contours;
}

// When computing the extent, ignore infinite, null, and NaN values
function finite(x) {
  return isFinite(x) ? x : NaN;
}

// Is the value above or equal to threshold?
// Handle null, undefined, and NaN as below threshold
function above(x, value) {
  return x != null && !isNaN(x = +x) ? x >= value : false;
}

// During smoothing, treat any invalid value as negative infinity
function valid(v) {
  return v == null || isNaN(v = +v) ? -Infinity : v;
}

// Linear interpolation for smoothing
function smooth1(x, v0, v1, value) {
  const a = value - v0;
  const b = v1 - v0;
  const d = isFinite(a) || isFinite(b) ? a / b : Math.sign(a) / Math.sign(b);
  return isNaN(d) ? x : x + d - 0.5;
}

function contourBand() {
  var dx = 1,
      dy = 1,
      threshold = d3Array.thresholdSturges,
      contourInstance = contours();
  
  function contourBand(values) {
    // Filter out null values for threshold calculation
    const validValues = values.filter(v => v != null && !isNaN(+v));
    if (validValues.length === 0) return []; // No valid data points
    
    var tz = threshold(validValues);

    // Convert number of thresholds into uniform thresholds
    if (!Array.isArray(tz)) {
      const e = d3Array.extent(validValues, finite);
      if (e[0] === undefined || e[1] === undefined) {
        return []; // No valid data
      }
      tz = d3Array.ticks(...d3Array.nice(e[0], e[1], tz), tz);
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

/**
 * 增强版等值线与等值面生成器
 * 基于Marching Squares算法，参考openhome.cc文章实现
 * 特点：1. 支持处理null值  2. 更平滑的等值线  3. 更好的边界处理  4. 鞍点二义性处理
 */

// 定义16种情况下的等值线连接方式
const CONTOUR_CASES = [
  [], // Case 0: 所有角点都低于阈值
  [[[1.0, 1.5], [0.5, 1.0]]], // Case 1
  [[[1.5, 1.0], [1.0, 1.5]]], // Case 2
  [[[1.5, 1.0], [0.5, 1.0]]], // Case 3
  [[[1.0, 0.5], [1.5, 1.0]]], // Case 4
  [[[1.0, 1.5], [0.5, 1.0]], [[1.0, 0.5], [1.5, 1.0]]], // Case 5: 鞍点
  [[[1.0, 0.5], [1.0, 1.5]]], // Case 6
  [[[1.0, 0.5], [0.5, 1.0]]], // Case 7
  [[[0.5, 1.0], [1.0, 0.5]]], // Case 8
  [[[1.0, 1.5], [1.0, 0.5]]], // Case 9
  [[[0.5, 1.0], [1.0, 0.5]], [[1.5, 1.0], [1.0, 1.5]]], // Case 10: 鞍点
  [[[1.5, 1.0], [1.0, 0.5]]], // Case 11
  [[[0.5, 1.0], [1.5, 1.0]]], // Case 12
  [[[1.0, 1.5], [1.5, 1.0]]], // Case 13
  [[[0.5, 1.0], [1.0, 1.5]]], // Case 14
  [] // Case 15: 所有角点都高于阈值
];

/**
 * 根据阈值将数据映射为角点，并标记其是否低于阈值
 * @param {Array<Array<Number>>} values 二维数据数组
 * @param {Number} threshold 阈值
 * @returns {Array<Array<Object>>} 角点数据
 */
function mapToCorners(values, threshold) {
  const corners = [];
  for (let r = 0; r < values.length; r++) {
    const row = [];
    for (let c = 0; c < values[r].length; c++) {
      const val = values[r][c];
      row.push({
        position: [c, r],
        value: val,
        belowThreshold: val != null && val < threshold
      });
    }
    corners.push(row);
  }
  return corners;
}

/**
 * 从角点数组创建细胞数据
 * @param {Array<Array<Object>>} corners 角点数组
 * @returns {Array<Object>} 细胞数组
 */
function createCells(corners) {
  const cells = [];
  for (let r = 0; r < corners.length - 1; r++) {
    for (let c = 0; c < corners[r].length - 1; c++) {
      // 提取单个细胞的四个角点
      const cellCorners = [
        corners[r][c],
        corners[r + 1][c],
        corners[r + 1][c + 1],
        corners[r][c + 1]
      ];
      
      // 计算细胞的边界框
      const bounds = {
        minX: c,
        minY: r,
        maxX: c + 1,
        maxY: r + 1
      };
      
      // 计算情况编码
      const caseCode = computeCaseCode(cellCorners);
      
      cells.push({
        corners: cellCorners,
        bounds: bounds,
        caseCode: caseCode
      });
    }
  }
  return cells;
}

/**
 * 计算细胞角点的情况编码
 * @param {Array<Object>} corners 细胞的四个角点
 * @returns {Number} 编码值(0-15)
 */
function computeCaseCode(corners) {
  return (corners[0].belowThreshold ? 8 : 0) +
         (corners[1].belowThreshold ? 1 : 0) + 
         (corners[2].belowThreshold ? 2 : 0) +
         (corners[3].belowThreshold ? 4 : 0);
}

/**
 * 在两点之间插值计算等值线交点的位置
 * @param {Object} p1 第一个点
 * @param {Object} p2 第二个点
 * @param {Number} threshold 阈值
 * @returns {Array} 插值后的位置
 */
function interpolatePoint(p1, p2, threshold) {
  // 处理null值的情况
  if (p1.value == null) return p2.position;
  if (p2.value == null) return p1.position;
  
  // 如果两个值相等，避免除以零的问题
  if (p1.value === p2.value) {
    return [(p1.position[0] + p2.position[0]) / 2, (p1.position[1] + p2.position[1]) / 2];
  }
  
  // 计算插值比例
  const t = (threshold - p1.value) / (p2.value - p1.value);
  
  // 应用平滑曲线插值，使曲线更加平滑（缓入缓出效果）
  const smoothT = t * t * (3 - 2 * t);
  
  // 线性插值计算位置
  return [
    p1.position[0] + smoothT * (p2.position[0] - p1.position[0]),
    p1.position[1] + smoothT * (p2.position[1] - p1.position[1])
  ];
}

/**
 * 计算四个角点的中心值
 * @param {Array<Object>} corners 四个角点
 * @returns {Number} 中心值
 */
function computeCenterValue(corners) {
  const validCorners = corners.filter(c => c.value != null && !isNaN(c.value));
  if (validCorners.length === 0) return -Infinity;
  
  // 计算平均值
  return validCorners.reduce((sum, c) => sum + c.value, 0) / validCorners.length;
}

/**
 * 生成细胞的等值线线段
 * @param {Object} cell 细胞数据
 * @param {Number} threshold 阈值
 * @returns {Array<Array<Array<Number>>>} 等值线线段数组
 */
function generateIsolines(cell, threshold) {
  const corners = cell.corners;
  const caseCode = cell.caseCode;
  
  // 处理鞍点情况（case 5和case 10），通过中心值来决定连接方式
  if (caseCode === 5 || caseCode === 10) {
    const centerValue = computeCenterValue(corners);
    if ((caseCode === 5 && centerValue >= threshold) || 
        (caseCode === 10 && centerValue >= threshold)) {
      // 使用替代连接方式
      return generateAlternativeSaddleIsolines(cell, threshold);
    }
  }
  
  // 获取此情况对应的线段模板
  const caseLines = CONTOUR_CASES[caseCode];
  if (!caseLines || caseLines.length === 0) return [];
  
  // 转换模板线段为实际坐标
  return caseLines.map(line => {
    return line.map(templatePoint => {
      const x = templatePoint[0];
      const y = templatePoint[1];
      
      // 确定要插值的两个角点
      let p1, p2;
      
      if (x === 0.5 && y === 1.0) { // 左边
        p1 = corners[0];
        p2 = corners[1];
      } else if (x === 1.0 && y === 1.5) { // 下边
        p1 = corners[1];
        p2 = corners[2];
      } else if (x === 1.5 && y === 1.0) { // 右边
        p1 = corners[2];
        p2 = corners[3];
      } else if (x === 1.0 && y === 0.5) { // 上边
        p1 = corners[3];
        p2 = corners[0];
      } else {
        // 不应该到这里，返回单元格中心点
        const center = [
          (corners[0].position[0] + corners[2].position[0]) / 2,
          (corners[0].position[1] + corners[2].position[1]) / 2
        ];
        return center;
      }
      
      return interpolatePoint(p1, p2, threshold);
    });
  });
}

/**
 * 生成鞍点情况的替代连接方式
 * @param {Object} cell 细胞数据
 * @param {Number} threshold 阈值
 * @returns {Array<Array<Array<Number>>>} 等值线线段数组
 */
function generateAlternativeSaddleIsolines(cell, threshold) {
  const corners = cell.corners;
  const caseCode = cell.caseCode;
  
  // 根据鞍点情况选择替代连接方式
  let alternativeCase;
  
  if (caseCode === 5) {
    // 替代情况: 连接左上-右下
    alternativeCase = [
      [[0.5, 1.0], [1.5, 1.0]]
    ];
  } else if (caseCode === 10) {
    // 替代情况: 连接左下-右上
    alternativeCase = [
      [[1.0, 0.5], [1.0, 1.5]]
    ];
  } else {
    return [];
  }
  
  // 转换替代连接方式为实际坐标
  return alternativeCase.map(line => {
    return line.map(templatePoint => {
      const x = templatePoint[0];
      const y = templatePoint[1];
      
      // 确定要插值的两个角点
      let p1, p2;
      
      if (x === 0.5 && y === 1.0) { // 左边
        p1 = corners[0];
        p2 = corners[1];
      } else if (x === 1.0 && y === 1.5) { // 下边
        p1 = corners[1];
        p2 = corners[2];
      } else if (x === 1.5 && y === 1.0) { // 右边
        p1 = corners[2];
        p2 = corners[3];
      } else if (x === 1.0 && y === 0.5) { // 上边
        p1 = corners[3];
        p2 = corners[0];
      } else {
        // 不应该到这里，返回单元格中心点
        const center = [
          (corners[0].position[0] + corners[2].position[0]) / 2,
          (corners[0].position[1] + corners[2].position[1]) / 2
        ];
        return center;
      }
      
      return interpolatePoint(p1, p2, threshold);
    });
  });
}

/**
 * 将线段拼接为闭合轮廓
 * @param {Array<Array<Array<Array<Number>>>>} lineSegments 线段数组
 * @returns {Array<Array<Array<Number>>>} 闭合轮廓数组
 */
function stitchContours(lineSegments) {
  // 扁平化所有线段
  const segments = lineSegments.flat();
  if (segments.length === 0) return [];
  
  // 存储已连接的轮廓和待处理的线段
  const contours = [];
  const remainingSegments = [...segments];
  
  // 当还有线段未处理时继续循环
  while (remainingSegments.length > 0) {
    // 取出一个线段作为起始线段
    const startSegment = remainingSegments.pop();
    const contour = [...startSegment];
    
    // 连接线段直到形成闭合轮廓或无法继续连接
    let connected = true;
    while (connected) {
      connected = false;
      
      // 获取当前轮廓的首尾点
      const firstPoint = contour[0];
      const lastPoint = contour[contour.length - 1];
      
      // 检查是否可以闭合轮廓
      if (pointsAreClose(firstPoint, lastPoint)) {
        // 轮廓已闭合，无需继续连接
        break;
      }
      
      // 尝试连接其他线段
      for (let i = 0; i < remainingSegments.length; i++) {
        const segment = remainingSegments[i];
        const segmentFirst = segment[0];
        const segmentLast = segment[segment.length - 1];
        
        // 检查是否可以连接到轮廓末尾
        if (pointsAreClose(lastPoint, segmentFirst)) {
          // 将线段除首点外的点添加到轮廓末尾
          for (let j = 1; j < segment.length; j++) {
            contour.push(segment[j]);
          }
          remainingSegments.splice(i, 1);
          connected = true;
          break;
        } 
        // 检查是否可以连接到轮廓开头
        else if (pointsAreClose(firstPoint, segmentLast)) {
          // 将线段除末点外的点添加到轮廓开头
          for (let j = segment.length - 2; j >= 0; j--) {
            contour.unshift(segment[j]);
          }
          remainingSegments.splice(i, 1);
          connected = true;
          break;
        }
        // 检查是否需要反转线段再连接
        else if (pointsAreClose(lastPoint, segmentLast)) {
          // 将反转的线段除首点外的点添加到轮廓末尾
          for (let j = segment.length - 2; j >= 0; j--) {
            contour.push(segment[j]);
          }
          remainingSegments.splice(i, 1);
          connected = true;
          break;
        }
        else if (pointsAreClose(firstPoint, segmentFirst)) {
          // 将反转的线段除末点外的点添加到轮廓开头
          for (let j = 1; j < segment.length; j++) {
            contour.unshift(segment[j]);
          }
          remainingSegments.splice(i, 1);
          connected = true;
          break;
        }
      }
    }
    
    // 应用平滑算法
    const smoothedContour = smoothContour(contour);
    
    // 添加到轮廓列表
    contours.push(smoothedContour);
  }
  
  return contours;
}

/**
 * 检查两个点是否接近（用于连接轮廓）
 * @param {Array<Number>} p1 第一个点
 * @param {Array<Number>} p2 第二个点
 * @returns {Boolean} 是否接近
 */
function pointsAreClose(p1, p2) {
  const EPSILON = 1e-10;
  const dx = p1[0] - p2[0];
  const dy = p1[1] - p2[1];
  return dx * dx + dy * dy < EPSILON;
}

/**
 * 平滑等值线轮廓
 * @param {Array<Array<Number>>} contour 轮廓点数组
 * @returns {Array<Array<Number>>} 平滑后的轮廓
 */
function smoothContour(contour) {
  if (contour.length <= 2) return contour;
  
  const result = [];
  const n = contour.length;
  
  // 拷贝首点
  result.push([...contour[0]]);
  
  // 对中间点应用Chaikin平滑算法
  for (let i = 0; i < n - 1; i++) {
    const p0 = contour[i];
    const p1 = contour[(i + 1) % n];
    
    // 避免在轮廓闭合处添加额外的点
    if (i === n - 2 && pointsAreClose(p1, contour[0])) {
      continue;
    }
    
    // 在两点之间插入平滑点
    const q = [
      0.75 * p0[0] + 0.25 * p1[0],
      0.75 * p0[1] + 0.25 * p1[1]
    ];
    const r = [
      0.25 * p0[0] + 0.75 * p1[0],
      0.25 * p0[1] + 0.75 * p1[1]
    ];
    
    result.push(q);
    result.push(r);
  }
  
  // 确保轮廓闭合
  if (!pointsAreClose(contour[n - 1], contour[0])) {
    result.push([...contour[n - 1]]);
  } else {
    result.push([...contour[0]]);
  }
  
  return result;
}

/**
 * 确定一个点是否在多边形内部
 * 这个函数用于getPointInsideRing函数的结果验证
 * @param {Array<Number>} point 要测试的点
 * @param {Array<Array<Number>>} polygon 多边形顶点数组
 * @returns {Boolean} 是否在内部
 */
function pointInPolygon(point, polygon) {
  let inside = false;
  const x = point[0], y = point[1];
  
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const xi = polygon[i][0], yi = polygon[i][1];
    const xj = polygon[j][0], yj = polygon[j][1];
    
    const intersect = ((yi > y) !== (yj > y)) && 
                      (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
    
    if (intersect) inside = !inside;
  }
  
  return inside;
}

/**
 * 获取一个环内部的点
 * @param {Array<Array<Number>>} ring 环的顶点数组
 * @returns {Array<Number>} 环内部的点
 */
function getPointInsideRing(ring) {
  if (ring.length < 3) {
    return null; // 少于3个点无法形成多边形
  }
  
  // 计算多边形的边界框
  let minX = Infinity, minY = Infinity;
  let maxX = -Infinity, maxY = -Infinity;
  
  for (const point of ring) {
    minX = Math.min(minX, point[0]);
    minY = Math.min(minY, point[1]);
    maxX = Math.max(maxX, point[0]);
    maxY = Math.max(maxY, point[1]);
  }
  
  // 计算中心点
  const centerX = (minX + maxX) / 2;
  const centerY = (minY + maxY) / 2;
  const centerPoint = [centerX, centerY];
  
  // 检查中心点是否在多边形内
  if (pointInPolygon(centerPoint, ring)) {
    return centerPoint;
  }
  
  // 如果中心点不在多边形内，尝试在边界框内随机生成点
  const MAX_ATTEMPTS = 50;
  for (let i = 0; i < MAX_ATTEMPTS; i++) {
    // 在边界框内随机生成点
    const randX = minX + Math.random() * (maxX - minX);
    const randY = minY + Math.random() * (maxY - minY);
    const testPoint = [randX, randY];
    
    if (pointInPolygon(testPoint, ring)) {
      return testPoint;
    }
  }
  
  // 尝试使用多边形边的中点并向内偏移
  for (let i = 0; i < ring.length - 1; i++) {
    const p1 = ring[i];
    const p2 = ring[(i + 1) % ring.length];
    
    // 计算边的中点
    const midX = (p1[0] + p2[0]) / 2;
    const midY = (p1[1] + p2[1]) / 2;
    
    // 计算边的法向量（向内）
    const dx = p2[0] - p1[0];
    const dy = p2[1] - p1[1];
    const length = Math.sqrt(dx * dx + dy * dy);
    
    if (length < 1e-10) continue; // 避免除以零
    
    // 法向量（逆时针旋转90度）
    const nx = -dy / length;
    const ny = dx / length;
    
    // 沿法向量向内偏移一小段距离
    const OFFSET = 0.1; // 小的偏移量
    const testPoint = [midX + nx * OFFSET, midY + ny * OFFSET];
    
    if (pointInPolygon(testPoint, ring)) {
      return testPoint;
    }
    
    // 尝试反方向
    const testPoint2 = [midX - nx * OFFSET, midY - ny * OFFSET];
    if (pointInPolygon(testPoint2, ring)) {
      return testPoint2;
    }
  }
  
  // 如果所有方法都失败，返回原始中心点（可能不准确）
  return centerPoint;
}

/**
 * 主等值线生成函数
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Array<Number>|Number} thresholds 阈值数组或单一阈值
 * @returns {Array<Object>|Object} GeoJSON格式的等值线数组或单个等值线
 */
function generateContours(data, thresholds) {
  if (!data || data.length === 0 || data[0].length === 0) {
    return { type: "MultiLineString", coordinates: [] };
  }
  
  // 如果传入的是单一阈值，转换为数组
  const isArray = Array.isArray(thresholds);
  const thresholdArray = isArray ? thresholds : [thresholds];
  
  // 为每个阈值生成等值线
  const results = thresholdArray.map(threshold => {
    // 1. 将数据映射为角点
    const corners = mapToCorners(data, threshold);
    
    // 2. 创建细胞数据
    const cells = createCells(corners);
    
    // 3. 为每个细胞生成等值线线段
    const lineSegments = cells.map(cell => generateIsolines(cell, threshold));
    
    // 4. 拼接线段为闭合轮廓
    const contours = stitchContours(lineSegments);
    
    return {
      type: "MultiLineString",
      threshold: threshold,
      coordinates: contours
    };
  });
  
  // 如果原始输入是单一阈值，返回单个结果，否则返回结果数组
  return isArray ? results : results[0];
}

/**
 * 生成等值面（带状区域）
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Array<Number>|Number} thresholds 阈值数组或单一阈值
 * @returns {Array<Object>|Object} GeoJSON格式的等值面数组或单个等值面
 */
function generateContourBands(data, thresholds) {
  if (!data || data.length === 0 || data[0].length === 0) {
    return { type: "MultiPolygon", coordinates: [] };
  }
  
  // 如果传入的是单一阈值，转换为数组
  const isArray = Array.isArray(thresholds);
  const thresholdArray = isArray ? thresholds : [thresholds];
  
  // 确保阈值数组已排序
  const sortedThresholds = [...thresholdArray].sort((a, b) => a - b);
  
  // 为每个阈值生成等值面
  const results = [];
  
  // 创建数据边界轮廓
  const width = data[0].length;
  const height = data.length;
  const boundaryRing = [
    [0, 0],
    [width - 1, 0],
    [width - 1, height - 1],
    [0, height - 1],
    [0, 0]
  ];
  
  // 获取所有等值线
  const allContours = generateContours(data, sortedThresholds);
  
  // 生成每个阈值的等值面
  for (let i = 0; i < sortedThresholds.length; i++) {
    const threshold = sortedThresholds[i];
    const contour = allContours[i];
    
    // 创建等值面
    const polygons = [];
    
    // 创建一个映射，标记数据点是否在阈值以上
    const isAboveThreshold = [];
    for (let y = 0; y < height; y++) {
      const row = [];
      for (let x = 0; x < width; x++) {
        const value = data[y][x];
        row.push(value != null && value >= threshold);
      }
      isAboveThreshold.push(row);
    }
    
    // 检查整个区域是否都在阈值的一侧
    let allAbove = true;
    let allBelow = true;
    
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        const value = data[y][x];
        if (value != null) {
          if (value < threshold) {
            allAbove = false;
          } else {
            allBelow = false;
          }
        }
        if (!allAbove && !allBelow) break;
      }
      if (!allAbove && !allBelow) break;
    }
    
    // 如果所有点都在阈值以上，创建一个覆盖整个区域的多边形
    if (allAbove) {
      polygons.push([boundaryRing]);
    }
    // 如果所有点都在阈值以下，不创建多边形（返回空）
    else if (!allBelow) {
      // 处理等值线轮廓
      if (contour.coordinates.length > 0) {
        for (const ring of contour.coordinates) {
          // 确保环是闭合的
          const closedRing = [...ring];
          if (closedRing.length > 0 && !pointsAreClose(closedRing[0], closedRing[closedRing.length - 1])) {
            closedRing.push([...closedRing[0]]);
          }
          
          // 确定这个轮廓是表示高于阈值的区域还是低于阈值的区域
          const testPoint = getPointInsideRing(closedRing);
          if (testPoint === null) continue; // 无法确定内部点，跳过
          
          const gridX = Math.floor(testPoint[0]);
          const gridY = Math.floor(testPoint[1]);
          
          // 检查这个点是否在网格范围内
          let isHighSide = false;
          if (gridX >= 0 && gridX < width && gridY >= 0 && gridY < height) {
            isHighSide = isAboveThreshold[gridY][gridX];
          }
          
          // 只保留高于阈值的区域
          if (isHighSide) {
            polygons.push([closedRing]);
          }
        }
      }
    }
    
    results.push({
      type: "MultiPolygon",
      threshold: threshold,
      coordinates: polygons
    });
  }
  
  // 如果原始输入是单一阈值，返回单个结果，否则返回结果数组
  return isArray ? results : results[0];
}

exports.contourBand = contourBand;
exports.contours = contours;
exports.generateContourBands = generateContourBands;
exports.generateContours = generateContours;

}));
