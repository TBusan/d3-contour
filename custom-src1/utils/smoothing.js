/**
 * 平滑处理工具
 * 用于平滑等值线
 */

/**
 * 平滑等值线
 * @param {Array<Array<Array<Number>>>} contours 等值线数组
 * @param {Number} smoothFactor 平滑因子（0-1）
 * @returns {Array<Array<Array<Number>>>} 平滑后的等值线
 */
export function smoothContours(contours, smoothFactor = 0.25) {
  if (!contours || contours.length === 0) {
    return contours;
  }
  
  // 确保平滑因子在有效范围内
  const factor = Math.max(0, Math.min(1, smoothFactor));
  
  // 对每个轮廓进行平滑
  return contours.map(contour => smoothContour(contour, factor));
}

/**
 * 平滑单个轮廓
 * @param {Array<Array<Number>>} contour 轮廓
 * @param {Number} factor 平滑因子
 * @returns {Array<Array<Number>>} 平滑后的轮廓
 */
function smoothContour(contour, factor) {
  if (contour.length < 3) {
    return contour;
  }
  
  // 检查轮廓是否闭合
  const isClosed = isContourClosed(contour);
  
  // 创建平滑后的轮廓
  const smoothed = [];
  
  // 对每个点进行平滑
  for (let i = 0; i < contour.length; i++) {
    const prev = getPrevPoint(contour, i, isClosed);
    const curr = contour[i];
    const next = getNextPoint(contour, i, isClosed);
    
    // 使用Chaikin平滑算法
    if (i === 0 && !isClosed) {
      // 如果是开放轮廓的第一个点，保持不变
      smoothed.push(curr);
    } else if (i === contour.length - 1 && !isClosed) {
      // 如果是开放轮廓的最后一个点，保持不变
      smoothed.push(curr);
    } else {
      // 计算当前点与前后点的中点
      const prevMid = [
        curr[0] * (1 - factor) + prev[0] * factor,
        curr[1] * (1 - factor) + prev[1] * factor
      ];
      
      const nextMid = [
        curr[0] * (1 - factor) + next[0] * factor,
        curr[1] * (1 - factor) + next[1] * factor
      ];
      
      // 添加两个中点
      smoothed.push(prevMid);
      smoothed.push(nextMid);
    }
  }
  
  // 如果是闭合轮廓，确保首尾相连
  if (isClosed) {
    smoothed.push(smoothed[0]);
  }
  
  return smoothed;
}

/**
 * 使用贝塞尔曲线平滑轮廓
 * @param {Array<Array<Number>>} contour 轮廓
 * @param {Number} factor 平滑因子
 * @returns {Array<Array<Number>>} 平滑后的轮廓
 */
export function bezierSmoothContour(contour, factor = 0.25) {
  if (contour.length < 3) {
    return contour;
  }
  
  // 检查轮廓是否闭合
  const isClosed = isContourClosed(contour);
  
  // 创建平滑后的轮廓
  const smoothed = [];
  
  // 对每个点生成贝塞尔曲线控制点
  for (let i = 0; i < contour.length - (isClosed ? 0 : 1); i++) {
    const p0 = contour[i];
    const p1 = contour[(i + 1) % contour.length];
    
    // 添加当前点
    smoothed.push(p0);
    
    // 如果不是最后一个点，添加贝塞尔曲线
    if (i < contour.length - 1 || isClosed) {
      const prev = getPrevPoint(contour, i, isClosed);
      const next = getNextPoint(contour, (i + 1) % contour.length, isClosed);
      
      // 计算控制点
      const cp1 = [
        p0[0] + (p1[0] - prev[0]) * factor,
        p0[1] + (p1[1] - prev[1]) * factor
      ];
      
      const cp2 = [
        p1[0] - (next[0] - p0[0]) * factor,
        p1[1] - (next[1] - p0[1]) * factor
      ];
      
      // 生成贝塞尔曲线上的点
      const bezierPoints = generateBezierPoints(p0, cp1, cp2, p1, 5);
      
      // 添加贝塞尔曲线上的点（除了起点和终点）
      for (let j = 1; j < bezierPoints.length - 1; j++) {
        smoothed.push(bezierPoints[j]);
      }
    }
  }
  
  // 如果是闭合轮廓，确保首尾相连
  if (isClosed) {
    smoothed.push(smoothed[0]);
  }
  
  return smoothed;
}

/**
 * 生成贝塞尔曲线上的点
 * @param {Array<Number>} p0 起点
 * @param {Array<Number>} p1 控制点1
 * @param {Array<Number>} p2 控制点2
 * @param {Array<Number>} p3 终点
 * @param {Number} numPoints 点的数量
 * @returns {Array<Array<Number>>} 贝塞尔曲线上的点
 */
function generateBezierPoints(p0, p1, p2, p3, numPoints) {
  const points = [];
  
  for (let i = 0; i <= numPoints; i++) {
    const t = i / numPoints;
    const point = cubicBezier(p0, p1, p2, p3, t);
    points.push(point);
  }
  
  return points;
}

/**
 * 计算三次贝塞尔曲线上的点
 * @param {Array<Number>} p0 起点
 * @param {Array<Number>} p1 控制点1
 * @param {Array<Number>} p2 控制点2
 * @param {Array<Number>} p3 终点
 * @param {Number} t 参数（0-1）
 * @returns {Array<Number>} 贝塞尔曲线上的点
 */
function cubicBezier(p0, p1, p2, p3, t) {
  const u = 1 - t;
  const tt = t * t;
  const uu = u * u;
  const uuu = uu * u;
  const ttt = tt * t;
  
  // B(t) = (1-t)^3 * P0 + 3(1-t)^2 * t * P1 + 3(1-t) * t^2 * P2 + t^3 * P3
  const x = uuu * p0[0] + 3 * uu * t * p1[0] + 3 * u * tt * p2[0] + ttt * p3[0];
  const y = uuu * p0[1] + 3 * uu * t * p1[1] + 3 * u * tt * p2[1] + ttt * p3[1];
  
  return [x, y];
}

/**
 * 获取前一个点
 * @param {Array<Array<Number>>} contour 轮廓
 * @param {Number} index 当前索引
 * @param {Boolean} isClosed 是否闭合
 * @returns {Array<Number>} 前一个点
 */
function getPrevPoint(contour, index, isClosed) {
  if (index === 0) {
    return isClosed ? contour[contour.length - 2] : contour[0];
  }
  return contour[index - 1];
}

/**
 * 获取后一个点
 * @param {Array<Array<Number>>} contour 轮廓
 * @param {Number} index 当前索引
 * @param {Boolean} isClosed 是否闭合
 * @returns {Array<Number>} 后一个点
 */
function getNextPoint(contour, index, isClosed) {
  if (index === contour.length - 1) {
    return isClosed ? contour[1] : contour[index];
  }
  return contour[index + 1];
}

/**
 * 检查轮廓是否闭合
 * @param {Array<Array<Number>>} contour 轮廓
 * @returns {Boolean} 是否闭合
 */
function isContourClosed(contour) {
  if (contour.length < 3) {
    return false;
  }
  
  const firstPoint = contour[0];
  const lastPoint = contour[contour.length - 1];
  
  const EPSILON = 1e-10;
  const dx = firstPoint[0] - lastPoint[0];
  const dy = firstPoint[1] - lastPoint[1];
  
  return dx * dx + dy * dy < EPSILON;
} 