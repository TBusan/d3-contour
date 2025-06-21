/**
 * 多边形辅助函数
 * 提供多边形操作的基础功能
 */

/**
 * 判断点是否在多边形内
 * @param {Array<Number>} point 点坐标 [x, y]
 * @param {Array<Array<Number>>} polygon 多边形
 * @returns {Boolean} 是否在多边形内
 */
export function isPointInPolygon(point, polygon) {
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

/**
 * 获取轮廓内的一个点
 * @param {Array<Array<Number>>} ring 轮廓
 * @returns {Array<Number>} 内部点
 */
export function getPointInsideRing(ring) {
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
 * 检查两个点是否接近
 * @param {Array<Number>} p1 第一个点
 * @param {Array<Number>} p2 第二个点
 * @returns {Boolean} 是否接近
 */
export function pointsAreClose(p1, p2) {
  const EPSILON = 1e-10;
  const dx = p1[0] - p2[0];
  const dy = p1[1] - p2[1];
  return dx * dx + dy * dy < EPSILON;
}

/**
 * 判断轮廓是否为顺时针方向
 * @param {Array<Array<Number>>} contour 轮廓
 * @returns {Boolean} 是否为顺时针
 */
export function isClockwise(contour) {
  let sum = 0;
  
  for (let i = 0; i < contour.length - 1; i++) {
    const p1 = contour[i];
    const p2 = contour[i + 1];
    sum += (p2[0] - p1[0]) * (p2[1] + p1[1]);
  }
  
  return sum > 0;
}

/**
 * 计算多边形面积
 * @param {Array<Array<Number>>} polygon 多边形
 * @returns {Number} 面积
 */
export function calculatePolygonArea(polygon) {
  let area = 0;
  
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    area += polygon[i][0] * polygon[j][1];
    area -= polygon[j][0] * polygon[i][1];
  }
  
  return Math.abs(area / 2);
}

/**
 * 计算两个多边形的交集
 * @param {Array<Array<Number>>} polygon1 多边形1
 * @param {Array<Array<Number>>} polygon2 多边形2
 * @returns {Array<Array<Number>>} 交集多边形
 */
export function polygonIntersection(polygon1, polygon2) {
  // 这个函数需要使用复杂的多边形裁剪算法
  // 如Sutherland-Hodgman算法或Weiler-Atherton算法
  // 这里只提供一个简化版本
  
  // 检查多边形1的点是否在多边形2内
  const points1Inside2 = polygon1.filter(p => isPointInPolygon(p, polygon2));
  
  // 检查多边形2的点是否在多边形1内
  const points2Inside1 = polygon2.filter(p => isPointInPolygon(p, polygon1));
  
  // 查找边的交点
  const intersections = [];
  
  for (let i = 0; i < polygon1.length; i++) {
    const i2 = (i + 1) % polygon1.length;
    const line1 = [polygon1[i], polygon1[i2]];
    
    for (let j = 0; j < polygon2.length; j++) {
      const j2 = (j + 1) % polygon2.length;
      const line2 = [polygon2[j], polygon2[j2]];
      
      const intersection = lineIntersection(line1, line2);
      if (intersection) {
        intersections.push(intersection);
      }
    }
  }
  
  // 合并所有点
  const allPoints = [...points1Inside2, ...points2Inside1, ...intersections];
  
  // 如果没有交点，返回空数组
  if (allPoints.length === 0) {
    return [];
  }
  
  // 按照凸包算法排序点
  return convexHull(allPoints);
}

/**
 * 计算两条线段的交点
 * @param {Array<Array<Number>>} line1 线段1
 * @param {Array<Array<Number>>} line2 线段2
 * @returns {Array<Number>|null} 交点或null
 */
function lineIntersection(line1, line2) {
  const x1 = line1[0][0];
  const y1 = line1[0][1];
  const x2 = line1[1][0];
  const y2 = line1[1][1];
  
  const x3 = line2[0][0];
  const y3 = line2[0][1];
  const x4 = line2[1][0];
  const y4 = line2[1][1];
  
  // 计算分母
  const denominator = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
  
  // 如果分母为0，线段平行或共线
  if (denominator === 0) {
    return null;
  }
  
  // 计算参数
  const ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denominator;
  const ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denominator;
  
  // 检查参数是否在[0,1]范围内（交点在两条线段上）
  if (ua >= 0 && ua <= 1 && ub >= 0 && ub <= 1) {
    // 计算交点坐标
    const x = x1 + ua * (x2 - x1);
    const y = y1 + ua * (y2 - y1);
    
    return [x, y];
  }
  
  return null;
}

/**
 * 计算点集的凸包
 * @param {Array<Array<Number>>} points 点集
 * @returns {Array<Array<Number>>} 凸包
 */
function convexHull(points) {
  // 如果点数小于3，无法形成凸包
  if (points.length < 3) {
    return points;
  }
  
  // 找到最左下角的点
  let startPoint = points[0];
  for (const point of points) {
    if (point[1] < startPoint[1] || (point[1] === startPoint[1] && point[0] < startPoint[0])) {
      startPoint = point;
    }
  }
  
  // 按照极角排序
  const sortedPoints = [...points].sort((a, b) => {
    if (a === startPoint) return -1;
    if (b === startPoint) return 1;
    
    const angleA = Math.atan2(a[1] - startPoint[1], a[0] - startPoint[0]);
    const angleB = Math.atan2(b[1] - startPoint[1], b[0] - startPoint[0]);
    
    if (angleA === angleB) {
      // 如果角度相同，选择距离更远的点
      const distA = (a[0] - startPoint[0]) ** 2 + (a[1] - startPoint[1]) ** 2;
      const distB = (b[0] - startPoint[0]) ** 2 + (b[1] - startPoint[1]) ** 2;
      return distA - distB;
    }
    
    return angleA - angleB;
  });
  
  // Graham扫描算法
  const hull = [sortedPoints[0], sortedPoints[1]];
  
  for (let i = 2; i < sortedPoints.length; i++) {
    while (hull.length > 1 && !isLeftTurn(hull[hull.length - 2], hull[hull.length - 1], sortedPoints[i])) {
      hull.pop();
    }
    hull.push(sortedPoints[i]);
  }
  
  return hull;
}

/**
 * 判断三个点是否形成左转
 * @param {Array<Number>} p1 点1
 * @param {Array<Number>} p2 点2
 * @param {Array<Number>} p3 点3
 * @returns {Boolean} 是否左转
 */
function isLeftTurn(p1, p2, p3) {
  return (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]) > 0;
} 