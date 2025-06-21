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
 * 为鞍点生成替代的等值线连接
 * @param {Object} cell 细胞数据
 * @param {Number} threshold 阈值
 * @returns {Array<Array<Array<Number>>>} 替代等值线线段
 */
function generateAlternativeSaddleIsolines(cell, threshold) {
  const corners = cell.corners;
  const caseCode = cell.caseCode;
  
  // 获取替代连接方式 (15 - 当前情况码)
  const altCaseLines = CONTOUR_CASES[15 - caseCode];
  if (!altCaseLines || altCaseLines.length === 0) return [];
  
  // 转换模板线段为实际坐标
  return altCaseLines.map(line => {
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
        // 中心点与边的连接，使用单元格中心
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
 * 合并线段片段为闭合轮廓
 * @param {Array<Array<Array<Number>>>} segments 线段片段
 * @returns {Array<Array<Array<Number>>>} 闭合轮廓
 */
function stitchContours(segments) {
  if (!segments || segments.length === 0) return [];
  
  // 展平所有线段
  let lineSegments = [];
  segments.forEach(cellSegments => {
    if (cellSegments) {
      cellSegments.forEach(segment => {
        if (segment && segment.length === 2) {
          lineSegments.push(segment);
        }
      });
    }
  });
  
  // 如果没有线段，直接返回
  if (lineSegments.length === 0) return [];
  
  // 构建轮廓线
  const contours = [];
  let currentContour = [lineSegments[0][0], lineSegments[0][1]];
  lineSegments.splice(0, 1);
  
  // 尝试闭合轮廓
  while (lineSegments.length > 0) {
    let foundMatch = false;
    let lastPoint = currentContour[currentContour.length - 1];
    
    // 查找能连接的下一个线段
    for (let i = 0; i < lineSegments.length; i++) {
      const segment = lineSegments[i];
      
      // 检查第一个点是否匹配
      if (pointsAreClose(lastPoint, segment[0])) {
        currentContour.push(segment[1]);
        lineSegments.splice(i, 1);
        foundMatch = true;
        break;
      }
      
      // 检查第二个点是否匹配
      if (pointsAreClose(lastPoint, segment[1])) {
        currentContour.push(segment[0]);
        lineSegments.splice(i, 1);
        foundMatch = true;
        break;
      }
    }
    
    // 如果无法继续连接，或者轮廓已闭合
    if (!foundMatch || pointsAreClose(currentContour[0], lastPoint)) {
      // 检查是否闭合（首尾相连）
      if (pointsAreClose(currentContour[0], lastPoint)) {
        // 确保首尾点完全一致
        currentContour[currentContour.length - 1] = [...currentContour[0]];
      }
      
      // 应用光滑处理
      const smoothedContour = smoothContour(currentContour);
      contours.push(smoothedContour);
      
      // 开始新的轮廓（如果还有剩余线段）
      if (lineSegments.length > 0) {
        currentContour = [lineSegments[0][0], lineSegments[0][1]];
        lineSegments.splice(0, 1);
      }
    }
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
 * 计算多边形区域的符号面积
 * @param {Array<Array<Number>>} ring 多边形顶点数组
 * @returns {Number} 符号面积
 */
function calculateArea(ring) {
  const n = ring.length;
  let area = 0;
  
  for (let i = 0; i < n; i++) {
    const j = (i + 1) % n;
    area += ring[i][0] * ring[j][1];
    area -= ring[j][0] * ring[i][1];
  }
  
  return area / 2;
}

/**
 * 确定一个点是否在多边形内部
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
 * 主等值线生成函数
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Number} threshold 阈值
 * @returns {Object} GeoJSON格式的等值线
 */
export function generateContours(data, threshold) {
  if (!data || data.length === 0 || data[0].length === 0) {
    return { type: "MultiLineString", coordinates: [] };
  }
  
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
}

/**
 * 生成等值面（带状区域）
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Number} lowerThreshold 下阈值
 * @param {Number} upperThreshold 上阈值
 * @returns {Object} GeoJSON格式的等值面
 */
export function generateContourBands(data, lowerThreshold, upperThreshold) {
  if (!data || data.length === 0 || data[0].length === 0) {
    return { type: "MultiPolygon", coordinates: [] };
  }
  
  // 获取下阈值的等值线
  const lowerContour = generateContours(data, lowerThreshold);
  
  // 获取上阈值的等值线
  const upperContour = generateContours(data, upperThreshold);
  
  // 构建多边形
  const polygons = [];
  
  // 处理下阈值轮廓作为外环
  lowerContour.coordinates.forEach(lowerRing => {
    // 确保环是闭合的
    if (!pointsAreClose(lowerRing[0], lowerRing[lowerRing.length - 1])) {
      lowerRing.push([...lowerRing[0]]);
    }
    
    // 计算面积以确保方向正确
    const area = calculateArea(lowerRing);
    
    // 如果面积为负，反转环的方向
    if (area < 0) {
      lowerRing.reverse();
    }
    
    // 查找此外环内的所有上阈值轮廓作为内环
    const holes = [];
    
    upperContour.coordinates.forEach(upperRing => {
      // 确保环是闭合的
      if (!pointsAreClose(upperRing[0], upperRing[upperRing.length - 1])) {
        upperRing.push([...upperRing[0]]);
      }
      
      // 检查上轮廓是否在下轮廓内部
      const testPoint = upperRing[0];
      if (pointInPolygon(testPoint, lowerRing)) {
        // 确保方向与外环相反
        const upperArea = calculateArea(upperRing);
        if (upperArea > 0) {
          upperRing.reverse();
        }
        
        holes.push(upperRing);
      }
    });
    
    // 构建多边形（带洞）
    polygons.push([lowerRing, ...holes]);
  });
  
  return {
    type: "MultiPolygon",
    lowerValue: lowerThreshold,
    upperValue: upperThreshold,
    coordinates: polygons
  };
} 