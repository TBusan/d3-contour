/**
 * Marching Squares算法实现
 * 用于生成等值线
 */

/**
 * Marching Squares算法生成等值线
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Number} threshold 阈值
 * @param {Object} options 配置选项
 * @returns {Array<Array<Array<Number>>>} 等值线数组
 */
export function marchingSquares(data, threshold, options = {}) {
  if (!data || data.length === 0 || data[0].length === 0) {
    return [];
  }
  
  const height = data.length;
  const width = data[0].length;
  
  // 存储生成的线段
  const segments = [];
  
  // 遍历每个单元格（2x2网格）
  for (let y = 0; y < height - 1; y++) {
    for (let x = 0; x < width - 1; x++) {
      // 获取单元格的四个角点值
      const topLeft = data[y][x];
      const topRight = data[y][x + 1];
      const bottomLeft = data[y + 1][x];
      const bottomRight = data[y + 1][x + 1];
      
      // 跳过包含null值的单元格
      if (topLeft == null || topRight == null || bottomLeft == null || bottomRight == null) {
        continue;
      }
      
      // 计算单元格类型（0-15）
      let cellType = 0;
      if (topLeft >= threshold) cellType |= 1;      // 0001
      if (topRight >= threshold) cellType |= 2;     // 0010
      if (bottomRight >= threshold) cellType |= 4;  // 0100
      if (bottomLeft >= threshold) cellType |= 8;   // 1000
      
      // 根据单元格类型生成线段
      const cellSegments = generateCellSegments(
        x, y, cellType, threshold,
        topLeft, topRight, bottomRight, bottomLeft,
        options
      );
      
      // 添加到线段列表
      segments.push(...cellSegments);
    }
  }
  
  // 连接线段形成完整的等值线
  return stitchSegments(segments, options.connectEnds !== false);
}

/**
 * 根据单元格类型生成线段
 * @param {Number} x 单元格x坐标
 * @param {Number} y 单元格y坐标
 * @param {Number} cellType 单元格类型（0-15）
 * @param {Number} threshold 阈值
 * @param {Number} topLeft 左上角值
 * @param {Number} topRight 右上角值
 * @param {Number} bottomRight 右下角值
 * @param {Number} bottomLeft 左下角值
 * @param {Object} options 配置选项
 * @returns {Array<Array<Array<Number>>>} 生成的线段
 */
function generateCellSegments(x, y, cellType, threshold, topLeft, topRight, bottomRight, bottomLeft, options) {
  // 根据Marching Squares查找表获取线段配置
  const config = MARCHING_SQUARES_LOOKUP[cellType];
  
  if (!config) {
    return [];
  }
  
  const segments = [];
  
  // 处理鞍点情况（类型5和10）
  if ((cellType === 5 || cellType === 10) && options && options.saddleResolution) {
    // 计算单元格中心值
    const centerValue = (topLeft + topRight + bottomRight + bottomLeft) / 4;
    
    // 根据中心值决定如何连接
    if (cellType === 5) {
      if (centerValue >= threshold) {
        // 连接左上-右下
        return generateCellSegments(x, y, 15, threshold, topLeft, topRight, bottomRight, bottomLeft, options);
      } else {
        // 连接右上-左下
        return generateCellSegments(x, y, 10, threshold, topLeft, topRight, bottomRight, bottomLeft, options);
      }
    } else { // cellType === 10
      if (centerValue >= threshold) {
        // 连接右上-左下
        return generateCellSegments(x, y, 15, threshold, topLeft, topRight, bottomRight, bottomLeft, options);
      } else {
        // 连接左上-右下
        return generateCellSegments(x, y, 5, threshold, topLeft, topRight, bottomRight, bottomLeft, options);
      }
    }
  }
  
  // 处理每个边的交点
  for (let i = 0; i < config.length; i += 2) {
    const edge1 = config[i];
    const edge2 = config[i + 1];
    
    // 计算第一个交点
    const point1 = calculateIntersection(
      x, y, edge1, threshold, topLeft, topRight, bottomRight, bottomLeft
    );
    
    // 计算第二个交点
    const point2 = calculateIntersection(
      x, y, edge2, threshold, topLeft, topRight, bottomRight, bottomLeft
    );
    
    // 添加线段
    segments.push([point1, point2]);
  }
  
  return segments;
}

/**
 * 计算等值线与单元格边的交点
 * @param {Number} x 单元格x坐标
 * @param {Number} y 单元格y坐标
 * @param {Number} edge 边的编号（0-3）
 * @param {Number} threshold 阈值
 * @param {Number} topLeft 左上角值
 * @param {Number} topRight 右上角值
 * @param {Number} bottomRight 右下角值
 * @param {Number} bottomLeft 左下角值
 * @returns {Array<Number>} 交点坐标 [x, y]
 */
function calculateIntersection(x, y, edge, threshold, topLeft, topRight, bottomRight, bottomLeft) {
  let t = 0;
  
  switch (edge) {
    case 0: // 上边
      t = (threshold - topLeft) / (topRight - topLeft);
      return [x + t, y];
    case 1: // 右边
      t = (threshold - topRight) / (bottomRight - topRight);
      return [x + 1, y + t];
    case 2: // 下边
      t = (threshold - bottomLeft) / (bottomRight - bottomLeft);
      return [x + t, y + 1];
    case 3: // 左边
      t = (threshold - topLeft) / (bottomLeft - topLeft);
      return [x, y + t];
    default:
      return [x, y]; // 不应该到这里
  }
}

/**
 * 将线段连接为完整的等值线
 * @param {Array<Array<Array<Number>>>} segments 线段数组
 * @param {Boolean} connectEnds 是否连接端点
 * @returns {Array<Array<Array<Number>>>} 连接后的等值线
 */
function stitchSegments(segments, connectEnds) {
  if (segments.length === 0) {
    return [];
  }
  
  // 创建一个副本，避免修改原始数据
  const remainingSegments = [...segments];
  const contours = [];
  
  // 连接线段直到没有剩余线段
  while (remainingSegments.length > 0) {
    // 取出一个线段作为起始线段
    const currentContour = [];
    let currentSegment = remainingSegments.pop();
    
    // 添加起始线段的点
    currentContour.push(currentSegment[0], currentSegment[1]);
    
    let connected = true;
    
    // 尝试连接更多线段
    while (connected && connectEnds) {
      connected = false;
      
      // 获取当前轮廓的首尾点
      const firstPoint = currentContour[0];
      const lastPoint = currentContour[currentContour.length - 1];
      
      // 尝试连接到其他线段
      for (let i = 0; i < remainingSegments.length; i++) {
        const segment = remainingSegments[i];
        const segmentStart = segment[0];
        const segmentEnd = segment[1];
        
        // 检查是否可以连接到轮廓末尾
        if (pointsAreClose(lastPoint, segmentStart)) {
          currentContour.push(segmentEnd);
          remainingSegments.splice(i, 1);
          connected = true;
          break;
        }
        // 检查是否可以连接到轮廓开头
        else if (pointsAreClose(firstPoint, segmentEnd)) {
          currentContour.unshift(segmentStart);
          remainingSegments.splice(i, 1);
          connected = true;
          break;
        }
        // 检查是否需要反转线段再连接
        else if (pointsAreClose(lastPoint, segmentEnd)) {
          currentContour.push(segmentStart);
          remainingSegments.splice(i, 1);
          connected = true;
          break;
        }
        else if (pointsAreClose(firstPoint, segmentStart)) {
          currentContour.unshift(segmentEnd);
          remainingSegments.splice(i, 1);
          connected = true;
          break;
        }
      }
      
      // 检查轮廓是否已闭合
      if (currentContour.length > 2 && pointsAreClose(currentContour[0], currentContour[currentContour.length - 1])) {
        connected = false; // 轮廓已闭合，停止连接
      }
    }
    
    // 添加到轮廓列表
    contours.push(currentContour);
  }
  
  return contours;
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
 * Marching Squares算法的查找表
 * 每个单元格类型对应的边交点配置
 * 格式: [edge1, edge2, edge3, edge4, ...]
 * 边的编号: 0=上边, 1=右边, 2=下边, 3=左边
 */
const MARCHING_SQUARES_LOOKUP = [
  [],             // 0: 0000 - 无交点
  [3, 0],         // 1: 0001 - 左上角
  [0, 1],         // 2: 0010 - 右上角
  [3, 1],         // 3: 0011 - 上边
  [1, 2],         // 4: 0100 - 右下角
  [3, 0, 1, 2],   // 5: 0101 - 左上角和右下角 (鞍点)
  [0, 2],         // 6: 0110 - 右边
  [3, 2],         // 7: 0111 - 左上角、右上角和右下角
  [2, 3],         // 8: 1000 - 左下角
  [2, 0],         // 9: 1001 - 左边
  [0, 1, 2, 3],   // 10: 1010 - 右上角和左下角 (鞍点)
  [2, 1],         // 11: 1011 - 左下角、左上角和右上角
  [1, 3],         // 12: 1100 - 下边
  [1, 0],         // 13: 1101 - 左下角、左上角和右下角
  [0, 3],         // 14: 1110 - 左下角、右上角和右下角
  []              // 15: 1111 - 无交点
]; 