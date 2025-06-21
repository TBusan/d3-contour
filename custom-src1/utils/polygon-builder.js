/**
 * 多边形生成工具
 * 用于将等值线转换为多边形
 */

/**
 * 从等值线创建多边形
 * @param {Array<Array<Array<Number>>>} contourLines 等值线数组
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Number} threshold 阈值
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Array<Array<Array<Number>>>} 多边形数组
 */
export function createPolygons(contourLines, data, threshold, width, height) {
  if (!contourLines || contourLines.length === 0) {
    // 如果没有等值线，创建边界多边形
    const boundaryPolygon = createDataBoundaryPolygon(width, height);
    
    // 确定边界多边形是否应该填充（基于阈值和数据）
    if (shouldFillBoundary(data, threshold)) {
      return [boundaryPolygon];
    } else {
      return []; // 不需要填充
    }
  }
  
  // 找到所有闭合的轮廓
  const closedContours = contourLines.filter(contour => {
    // 检查首尾点是否接近
    if (contour.length < 3) return false;
    
    const firstPoint = contour[0];
    const lastPoint = contour[contour.length - 1];
    
    return pointsAreClose(firstPoint, lastPoint);
  });
  
  // 如果没有闭合轮廓，尝试闭合它们
  if (closedContours.length === 0) {
    const boundaryContours = createBoundaryContours(contourLines, data, threshold, width, height);
    return boundaryContours;
  }
  
  // 构建多边形层次结构（确定哪些轮廓在其他轮廓内）
  const hierarchy = buildContourHierarchy(closedContours);
  
  // 将层次结构转换为GeoJSON多边形
  return convertHierarchyToPolygons(hierarchy);
}

/**
 * 检查两个点是否接近
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
 * 创建数据边界多边形
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Array<Array<Number>>} 边界多边形
 */
function createDataBoundaryPolygon(width, height) {
  // 创建边界轮廓
  const boundaryPoints = [];
  
  // 添加四个角点（按顺时针方向）
  boundaryPoints.push([0, 0]);
  boundaryPoints.push([width - 1, 0]);
  boundaryPoints.push([width - 1, height - 1]);
  boundaryPoints.push([0, height - 1]);
  boundaryPoints.push([0, 0]); // 闭合轮廓
  
  return boundaryPoints;
}

/**
 * 判断边界是否应该填充
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Number} threshold 阈值
 * @returns {Boolean} 是否应该填充
 */
function shouldFillBoundary(data, threshold) {
  // 检查角点和中心点的值
  const height = data.length;
  const width = data[0].length;
  
  // 检查四个角点
  const cornerValues = [
    data[0][0],
    data[0][width - 1],
    data[height - 1][0],
    data[height - 1][width - 1]
  ];
  
  // 检查中心点
  const centerX = Math.floor(width / 2);
  const centerY = Math.floor(height / 2);
  const centerValue = data[centerY][centerX];
  
  // 如果大多数边界点的值大于阈值，则填充边界
  const validValues = cornerValues.filter(v => v != null);
  if (validValues.length === 0) return false;
  
  const aboveThresholdCount = validValues.filter(v => v >= threshold).length;
  return aboveThresholdCount > validValues.length / 2 || (centerValue != null && centerValue >= threshold);
}

/**
 * 创建边界轮廓
 * @param {Array<Array<Array<Number>>>} contourLines 等值线数组
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Number} threshold 阈值
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Array<Array<Array<Number>>>} 多边形数组
 */
function createBoundaryContours(contourLines, data, threshold, width, height) {
  // 创建数据边界多边形
  const boundaryPoints = createDataBoundaryPolygon(width, height);
  
  // 将不闭合的等值线与边界连接起来
  const processedContours = [];
  const openContours = [];
  
  // 找到所有开放的轮廓（不闭合的）
  for (const contour of contourLines) {
    if (contour.length < 2) continue;
    
    const firstPoint = contour[0];
    const lastPoint = contour[contour.length - 1];
    
    if (pointsAreClose(firstPoint, lastPoint)) {
      // 闭合轮廓，直接添加
      processedContours.push([...contour]);
    } else {
      // 开放轮廓，需要与边界连接
      openContours.push([...contour]);
    }
  }
  
     // 处理开放轮廓
   if (openContours.length > 0) {
     // 将开放轮廓两两配对，形成闭合轮廓
     const pairedContours = pairOpenContours(openContours);
    
    for (const pairedContour of pairedContours) {
      // 创建一个新的闭合轮廓
      processedContours.push(pairedContour);
    }
    
    // 处理剩余的未配对轮廓
    for (const openContour of openContours) {
      if (openContour.processed) continue;
      
      // 将未配对的轮廓与边界连接
      const closedContour = connectContourToBoundary(openContour, width, height);
      if (closedContour.length > 2) {
        processedContours.push(closedContour);
      }
    }
  }
  
  // 确定哪些轮廓是"洞"
  const holes = [];
  const shells = [];
  
  for (const contour of processedContours) {
    // 确定这个轮廓是外壳还是洞
    const isHole = isContourHole(contour, data, threshold);
    
    // 确保轮廓是闭合的
    const closedContour = [...contour];
    if (!pointsAreClose(closedContour[0], closedContour[closedContour.length - 1])) {
      closedContour.push(closedContour[0]); // 闭合轮廓
    }
    
    if (isHole) {
      holes.push(closedContour);
    } else {
      shells.push(closedContour);
    }
  }
  
  // 创建多边形
  const polygons = [];
  
  // 如果没有外壳，检查是否应该使用边界作为外壳
  if (shells.length === 0 && shouldFillBoundary(data, threshold)) {
    const polygon = [boundaryPoints, ...holes];
    polygons.push(polygon);
  } else {
    // 为每个外壳创建一个多边形，并找到其中的洞
    for (const shell of shells) {
      const shellHoles = holes.filter(hole => isPointInPolygon(getPointInsideRing(hole), shell));
      const polygon = [shell, ...shellHoles];
      polygons.push(polygon);
    }
  }
  
  return polygons;
}

/**
 * 将开放轮廓两两配对
 * @param {Array<Array<Array<Number>>>} openContours 开放轮廓数组
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Array<Array<Array<Number>>>} 配对后的闭合轮廓
 */
function pairOpenContours(openContours) {
  const pairedContours = [];
  
  // 标记所有轮廓为未处理
  for (const contour of openContours) {
    contour.processed = false;
  }
  
  // 尝试配对轮廓
  for (let i = 0; i < openContours.length; i++) {
    if (openContours[i].processed) continue;
    
    const contour1 = openContours[i];
    const end1 = contour1[contour1.length - 1];
    
    let bestMatchIndex = -1;
    let minDistance = Infinity;
    
    // 找到最佳配对轮廓
    for (let j = 0; j < openContours.length; j++) {
      if (i === j || openContours[j].processed) continue;
      
      const contour2 = openContours[j];
      const start2 = contour2[0];
      const end2 = contour2[contour2.length - 1];
      
      // 计算端点之间的距离
      const d1 = distance(end1, start2);
      const d2 = distance(end1, end2);
      
      // 选择最小距离
      if (d1 < minDistance) {
        minDistance = d1;
        bestMatchIndex = j;
      }
      
      if (d2 < minDistance) {
        minDistance = d2;
        bestMatchIndex = j;
      }
    }
    
    // 如果找到合适的配对
    if (bestMatchIndex >= 0 && minDistance < 5.0) { // 设置合理的阈值
      const contour2 = openContours[bestMatchIndex];
      const start2 = contour2[0];
      const end2 = contour2[contour2.length - 1];
      
      // 创建一个新的闭合轮廓
      let pairedContour;
      
      if (distance(end1, start2) < distance(end1, end2)) {
        // end1连接到start2
        pairedContour = [...contour1, ...contour2];
      } else {
        // end1连接到end2，需要反转contour2
        pairedContour = [...contour1, ...contour2.slice().reverse()];
      }
      
      pairedContours.push(pairedContour);
      
      // 标记为已处理
      contour1.processed = true;
      contour2.processed = true;
    }
  }
  
  return pairedContours;
}

/**
 * 计算两点之间的距离
 * @param {Array<Number>} p1 第一个点
 * @param {Array<Number>} p2 第二个点
 * @returns {Number} 距离
 */
function distance(p1, p2) {
  const dx = p1[0] - p2[0];
  const dy = p1[1] - p2[1];
  return Math.sqrt(dx * dx + dy * dy);
}

/**
 * 将开放轮廓与数据边界连接
 * @param {Array<Array<Number>>} contour 开放轮廓
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Array<Array<Number>>} 闭合轮廓
 */
function connectContourToBoundary(contour, width, height) {
  if (contour.length < 2) return contour;
  
  const start = contour[0];
  const end = contour[contour.length - 1];
  
  // 确定起点和终点是否靠近边界
  const startNearBoundary = isPointNearBoundary(start, width, height);
  const endNearBoundary = isPointNearBoundary(end, width, height);
  
  if (!startNearBoundary && !endNearBoundary) {
    return contour; // 两端都不靠近边界，无法闭合
  }
  
  // 创建一个新的闭合轮廓
  const closedContour = [...contour];
  
  // 计算起点和终点到各边界的最短距离
  const startBoundary = findNearestBoundaryPoint(start, width, height);
  const endBoundary = findNearestBoundaryPoint(end, width, height);
  
  // 添加边界点以闭合轮廓
  if (startNearBoundary && endNearBoundary) {
    // 两端都靠近边界，沿着边界连接它们
    const boundaryPath = createBoundaryPath(endBoundary, startBoundary, width, height);
    closedContour.push(...boundaryPath);
  } else if (startNearBoundary) {
    // 只有起点靠近边界，尝试将终点也连接到边界
    closedContour.push(startBoundary);
  } else if (endNearBoundary) {
    // 只有终点靠近边界，尝试将起点也连接到边界
    closedContour.unshift(endBoundary);
  }
  
  return closedContour;
}

/**
 * 判断点是否靠近数据边界
 * @param {Array<Number>} point 点坐标
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Boolean} 是否靠近边界
 */
function isPointNearBoundary(point, width, height) {
  const [x, y] = point;
  const THRESHOLD = 0.01; // 边界阈值
  
  return x < THRESHOLD || x > width - 1 - THRESHOLD || 
         y < THRESHOLD || y > height - 1 - THRESHOLD;
}

/**
 * 找到最近的边界点
 * @param {Array<Number>} point 点坐标
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Array<Number>} 边界点坐标
 */
function findNearestBoundaryPoint(point, width, height) {
  const [x, y] = point;
  
  // 计算到各边界的距离
  const distToLeft = x;
  const distToRight = width - 1 - x;
  const distToTop = y;
  const distToBottom = height - 1 - y;
  
  // 找到最近的边界
  const minDist = Math.min(distToLeft, distToRight, distToTop, distToBottom);
  
  if (minDist === distToLeft) {
    return [0, y]; // 左边界
  } else if (minDist === distToRight) {
    return [width - 1, y]; // 右边界
  } else if (minDist === distToTop) {
    return [x, 0]; // 上边界
  } else {
    return [x, height - 1]; // 下边界
  }
}

/**
 * 创建沿边界的路径
 * @param {Array<Number>} start 起点
 * @param {Array<Number>} end 终点
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Array<Array<Number>>} 边界路径
 */
function createBoundaryPath(start, end, width, height) {
  const path = [];
  
  // 确定起点和终点在哪个边界上
  const startEdge = getBoundaryEdge(start, width, height);
  const endEdge = getBoundaryEdge(end, width, height);
  
  if (startEdge === endEdge) {
    // 如果在同一边界上，直接连接
    return [];
  }
  
  // 按顺时针方向添加角点
  const corners = [
    [0, 0],
    [width - 1, 0],
    [width - 1, height - 1],
    [0, height - 1]
  ];
  
  // 确定起点后的第一个角点索引
  let currentEdge = startEdge;
  
  // 按顺时针方向遍历边界，直到到达终点所在的边
  while (currentEdge !== endEdge) {
    const nextEdge = (currentEdge + 1) % 4;
    const cornerPoint = corners[nextEdge];
    
    path.push(cornerPoint);
    currentEdge = nextEdge;
  }
  
  return path;
}

/**
 * 获取点所在的边界
 * @param {Array<Number>} point 点坐标
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Number} 边界索引（0=上，1=右，2=下，3=左）
 */
function getBoundaryEdge(point, width, height) {
  const [x, y] = point;
  const EPSILON = 1e-10;
  
  if (Math.abs(y) < EPSILON) return 0; // 上边界
  if (Math.abs(x - (width - 1)) < EPSILON) return 1; // 右边界
  if (Math.abs(y - (height - 1)) < EPSILON) return 2; // 下边界
  if (Math.abs(x) < EPSILON) return 3; // 左边界
  
  // 如果不在边界上，返回最近的边界
  const distToTop = y;
  const distToRight = width - 1 - x;
  const distToBottom = height - 1 - y;
  const distToLeft = x;
  
  const minDist = Math.min(distToTop, distToRight, distToBottom, distToLeft);
  
  if (minDist === distToTop) return 0;
  if (minDist === distToRight) return 1;
  if (minDist === distToBottom) return 2;
  return 3;
}

/**
 * 构建轮廓层次结构
 * @param {Array<Array<Array<Number>>>} contours 轮廓数组
 * @returns {Array<Object>} 层次结构
 */
function buildContourHierarchy(contours) {
  const nodes = contours.map(contour => ({
    contour,
    children: []
  }));
  
  // 构建包含关系
  for (let i = 0; i < nodes.length; i++) {
    const nodeI = nodes[i];
    
    for (let j = 0; j < nodes.length; j++) {
      if (i === j) continue;
      
      const nodeJ = nodes[j];
      
      // 检查nodeJ是否包含在nodeI中
      const pointInJ = getPointInsideRing(nodeJ.contour);
      if (isPointInPolygon(pointInJ, nodeI.contour)) {
        // 检查是否已经有更近的父节点
        let hasCloserParent = false;
        
        for (let k = 0; k < nodes.length; k++) {
          if (k === i || k === j) continue;
          
          const nodeK = nodes[k];
          
          if (isPointInPolygon(pointInJ, nodeK.contour) && 
              isPointInPolygon(getPointInsideRing(nodeK.contour), nodeI.contour)) {
            hasCloserParent = true;
            break;
          }
        }
        
        if (!hasCloserParent) {
          nodeI.children.push(nodeJ);
        }
      }
    }
  }
  
  // 找到根节点（没有父节点的节点）
  const rootNodes = nodes.filter(node => {
    for (const otherNode of nodes) {
      if (node !== otherNode && 
          isPointInPolygon(getPointInsideRing(node.contour), otherNode.contour)) {
        return false;
      }
    }
    return true;
  });
  
  return rootNodes;
}

/**
 * 将层次结构转换为多边形
 * @param {Array<Object>} hierarchy 层次结构
 * @returns {Array<Array<Array<Number>>>} 多边形数组
 */
function convertHierarchyToPolygons(hierarchy) {
  const polygons = [];
  
  for (const node of hierarchy) {
    const shell = node.contour;
    const holes = node.children.map(child => child.contour);
    
    // 确保轮廓是闭合的
    const closedShell = [...shell];
    if (!pointsAreClose(closedShell[0], closedShell[closedShell.length - 1])) {
      closedShell.push(closedShell[0]);
    }
    
    const closedHoles = holes.map(hole => {
      const closedHole = [...hole];
      if (!pointsAreClose(closedHole[0], closedHole[closedHole.length - 1])) {
        closedHole.push(closedHole[0]);
      }
      return closedHole;
    });
    
    // 创建多边形
    const polygon = [closedShell, ...closedHoles];
    polygons.push(polygon);
    
    // 递归处理子节点的子节点
    for (const child of node.children) {
      const childPolygons = convertHierarchyToPolygons(child.children);
      polygons.push(...childPolygons);
    }
  }
  
  return polygons;
}

/**
 * 判断轮廓是否为"洞"
 * @param {Array<Array<Number>>} contour 轮廓
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Number} threshold 阈值
 * @returns {Boolean} 是否为洞
 */
function isContourHole(contour, data, threshold) {
  // 获取轮廓内的一个点
  const insidePoint = getPointInsideRing(contour);
  
  // 找到最近的数据点
  const x = Math.round(insidePoint[0]);
  const y = Math.round(insidePoint[1]);
  
  // 检查该点的值是否小于阈值
  if (x >= 0 && x < data[0].length && y >= 0 && y < data.length) {
    return data[y][x] < threshold;
  }
  
  // 如果无法确定，使用轮廓方向
  return isClockwise(contour);
}

/**
 * 判断轮廓是否为顺时针方向
 * @param {Array<Array<Number>>} contour 轮廓
 * @returns {Boolean} 是否为顺时针
 */
function isClockwise(contour) {
  let sum = 0;
  
  for (let i = 0; i < contour.length - 1; i++) {
    const p1 = contour[i];
    const p2 = contour[i + 1];
    sum += (p2[0] - p1[0]) * (p2[1] + p1[1]);
  }
  
  return sum > 0;
}

/**
 * 获取轮廓内的一个点
 * @param {Array<Array<Number>>} ring 轮廓
 * @returns {Array<Number>} 内部点
 */
function getPointInsideRing(ring) {
  // 计算轮廓的中心点
  let sumX = 0;
  let sumY = 0;
  
  for (const point of ring) {
    sumX += point[0];
    sumY += point[1];
  }
  
  const centerX = sumX / ring.length;
  const centerY = sumY / ring.length;
  
  // 检查中心点是否在轮廓内
  if (isPointInPolygon([centerX, centerY], ring)) {
    return [centerX, centerY];
  }
  
  // 如果中心点不在轮廓内，尝试找到一个在轮廓内的点
  // 使用光线投射法的变体
  for (let i = 0; i < ring.length - 1; i++) {
    const p1 = ring[i];
    const p2 = ring[i + 1];
    
    // 计算边的中点
    const midX = (p1[0] + p2[0]) / 2;
    const midY = (p1[1] + p2[1]) / 2;
    
    // 计算边的法向量（向内）
    const dx = p2[0] - p1[0];
    const dy = p2[1] - p1[1];
    const length = Math.sqrt(dx * dx + dy * dy);
    
    if (length > 0) {
      // 单位法向量
      const nx = -dy / length;
      const ny = dx / length;
      
      // 沿法向量移动一小段距离
      const insideX = midX + nx * 0.01;
      const insideY = midY + ny * 0.01;
      
      // 检查这个点是否在轮廓内
      if (isPointInPolygon([insideX, insideY], ring)) {
        return [insideX, insideY];
      }
      
      // 尝试相反方向
      const outsideX = midX - nx * 0.01;
      const outsideY = midY - ny * 0.01;
      
      if (isPointInPolygon([outsideX, outsideY], ring)) {
        return [outsideX, outsideY];
      }
    }
  }
  
  // 如果上述方法都失败，返回轮廓的第一个点（不理想但至少有一个点）
  return ring[0];
}

/**
 * 判断点是否在多边形内
 * @param {Array<Number>} point 点坐标 [x, y]
 * @param {Array<Array<Number>>} polygon 多边形
 * @returns {Boolean} 是否在多边形内
 */
function isPointInPolygon(point, polygon) {
  // 使用光线投射算法
  let inside = false;
  const x = point[0];
  const y = point[1];
  
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const xi = polygon[i][0];
    const yi = polygon[i][1];
    const xj = polygon[j][0];
    const yj = polygon[j][1];
    
    const intersect = ((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
    if (intersect) inside = !inside;
  }
  
  return inside;
} 