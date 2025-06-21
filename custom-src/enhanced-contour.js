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
  
  // 创建等值面
  const polygons = [];
  
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
  
  // 创建一个映射，标记数据点是否在下阈值以上
  const isAboveLower = [];
  for (let y = 0; y < height; y++) {
    const row = [];
    for (let x = 0; x < width; x++) {
      const value = data[y][x];
      row.push(value != null && value >= lowerThreshold);
    }
    isAboveLower.push(row);
  }
  
  // 创建一个映射，标记数据点是否在上阈值以下
  const isBelowUpper = [];
  for (let y = 0; y < height; y++) {
    const row = [];
    for (let x = 0; x < width; x++) {
      const value = data[y][x];
      row.push(value != null && value < upperThreshold);
    }
    isBelowUpper.push(row);
  }
  
  // 处理下阈值轮廓作为外环
  if (lowerContour.coordinates.length > 0) {
    for (const lowerRing of lowerContour.coordinates) {
      // 确保环是闭合的
      const closedLowerRing = [...lowerRing];
      if (closedLowerRing.length > 0 && !pointsAreClose(closedLowerRing[0], closedLowerRing[closedLowerRing.length - 1])) {
        closedLowerRing.push([...closedLowerRing[0]]);
      }
      
      // 计算面积以确保方向正确（外环应为顺时针）
      const area = calculateArea(closedLowerRing);
      if (area < 0) {
        // 如果是逆时针，反转为顺时针
        closedLowerRing.reverse();
      }
      
      // 创建多边形，初始只有外环
      const polygon = [closedLowerRing];
      
      // 查找此外环内的所有上阈值轮廓作为内环
      for (const upperRing of upperContour.coordinates) {
        // 确保环是闭合的
        const closedUpperRing = [...upperRing];
        if (closedUpperRing.length > 0 && !pointsAreClose(closedUpperRing[0], closedUpperRing[closedUpperRing.length - 1])) {
          closedUpperRing.push([...closedUpperRing[0]]);
        }
        
        // 计算面积以确保方向正确（内环应为逆时针）
        const upperArea = calculateArea(closedUpperRing);
        if (upperArea > 0) {
          // 如果是顺时针，反转为逆时针
          closedUpperRing.reverse();
        }
        
        // 检查上轮廓是否在下轮廓内部
        if (closedUpperRing.length > 0) {
          // 使用多个点来确定是否在内部，增加可靠性
          let insideCount = 0;
          const testPoints = [
            closedUpperRing[0],
            closedUpperRing[Math.floor(closedUpperRing.length / 3)],
            closedUpperRing[Math.floor(closedUpperRing.length * 2 / 3)]
          ];
          
          for (const testPoint of testPoints) {
            if (pointInPolygon(testPoint, closedLowerRing)) {
              insideCount++;
            }
          }
          
          // 如果大多数测试点在内部，则认为是内环
          if (insideCount >= 2) {
            polygon.push(closedUpperRing);
          }
        }
      }
      
      polygons.push(polygon);
    }
  }
  
  // 处理边界情况：如果没有下阈值轮廓，或者需要处理外部区域
  if (polygons.length === 0) {
    // 检查是否整个区域都在阈值范围内
    let allInRange = true;
    let anyInRange = false;
    
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        const value = data[y][x];
        if (value != null) {
          if (value >= lowerThreshold && value < upperThreshold) {
            anyInRange = true;
          } else {
            allInRange = false;
          }
        }
      }
    }
    
    // 如果整个区域都在范围内，使用边界作为外环
    if (allInRange || anyInRange) {
      // 创建边界多边形
      const boundaryPolygon = [boundaryRing];
      
      // 添加所有上阈值轮廓作为内环
      for (const upperRing of upperContour.coordinates) {
        // 确保环是闭合的
        const closedUpperRing = [...upperRing];
        if (closedUpperRing.length > 0 && !pointsAreClose(closedUpperRing[0], closedUpperRing[closedUpperRing.length - 1])) {
          closedUpperRing.push([...closedUpperRing[0]]);
        }
        
        // 确保内环为逆时针方向
        const upperArea = calculateArea(closedUpperRing);
        if (upperArea > 0) {
          closedUpperRing.reverse();
        }
        
        boundaryPolygon.push(closedUpperRing);
      }
      
      polygons.push(boundaryPolygon);
    }
  }
  
  // 如果仍然没有多边形，检查是否有特殊情况需要处理
  if (polygons.length === 0 && upperContour.coordinates.length > 0) {
    // 尝试使用上阈值轮廓的补集作为等值面
    for (const upperRing of upperContour.coordinates) {
      // 确保环是闭合的
      const closedUpperRing = [...upperRing];
      if (closedUpperRing.length > 0 && !pointsAreClose(closedUpperRing[0], closedUpperRing[closedUpperRing.length - 1])) {
        closedUpperRing.push([...closedUpperRing[0]]);
      }
      
      // 计算面积以确保方向正确
      const area = calculateArea(closedUpperRing);
      if (area > 0) {
        // 如果是顺时针，反转为逆时针（作为内环）
        closedUpperRing.reverse();
      }
      
      // 创建一个使用边界作为外环，上阈值轮廓作为内环的多边形
      polygons.push([boundaryRing, closedUpperRing]);
    }
  }
  
  return {
    type: "MultiPolygon",
    lowerValue: lowerThreshold,
    upperValue: upperThreshold,
    coordinates: polygons
  };
}

/**
 * 同时生成等值线和等值面
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Number} lowerThreshold 下阈值
 * @param {Number} upperThreshold 上阈值
 * @returns {Object} 包含等值线和等值面的对象
 */
export function generateContourAndBands(data, lowerThreshold, upperThreshold) {
  const contours = generateContours(data, lowerThreshold);
  const bands = generateContourBands(data, lowerThreshold, upperThreshold);
  
  return {
    contours,
    bands
  };
} 