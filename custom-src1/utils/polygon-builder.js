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
    return [];
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
 * 创建边界轮廓
 * @param {Array<Array<Array<Number>>>} contourLines 等值线数组
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Number} threshold 阈值
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Array<Array<Array<Number>>>} 多边形数组
 */
function createBoundaryContours(contourLines, data, threshold, width, height) {
  // 创建边界轮廓
  const boundaryPoints = [];
  
  // 添加四个角点
  boundaryPoints.push([0, 0]);
  boundaryPoints.push([width - 1, 0]);
  boundaryPoints.push([width - 1, height - 1]);
  boundaryPoints.push([0, height - 1]);
  boundaryPoints.push([0, 0]); // 闭合轮廓
  
  // 创建一个包含所有点的多边形
  const boundaryPolygon = [boundaryPoints];
  
  // 如果没有等值线，返回边界多边形
  if (contourLines.length === 0) {
    return [boundaryPolygon];
  }
  
  // 确定哪些轮廓是"洞"
  const holes = [];
  const shells = [];
  
  for (const contour of contourLines) {
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
  
  // 如果没有外壳，使用边界作为外壳
  if (shells.length === 0) {
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