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
export function generateContours(data, thresholds) {
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
export function generateContourBands(data, thresholds) {
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

/**
 * 同时生成等值线和等值面
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Array<Number>|Number} thresholds 阈值数组或单一阈值
 * @returns {Object} 包含等值线和等值面的对象
 */
export function generateContourAndBands(data, thresholds) {
  const contours = generateContours(data, thresholds);
  const bands = generateContourBands(data, thresholds);
  
  return {
    contours,
    bands
  };
} 