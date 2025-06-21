/**
 * 多边形工具类
 * 用于处理多边形操作
 */

import { isPointInPolygon, getPointInsideRing } from './polygon-helpers.js';

/**
 * 合并多边形
 * @param {Array<Array<Array<Number>>>} lowerPolygons 低阈值多边形
 * @param {Array<Array<Array<Number>>>} upperPolygons 高阈值多边形
 * @returns {Array<Array<Array<Number>>>} 合并后的多边形
 */
export function mergePolygons(lowerPolygons, upperPolygons) {
  // 如果没有低阈值多边形，直接返回高阈值多边形
  if (!lowerPolygons || lowerPolygons.length === 0) {
    return upperPolygons || [];
  }
  
  // 如果没有高阈值多边形，直接返回低阈值多边形
  if (!upperPolygons || upperPolygons.length === 0) {
    return lowerPolygons;
  }
  
  // 创建结果多边形
  const resultPolygons = [];
  
  // 处理每个低阈值多边形
  for (const lowerPolygon of lowerPolygons) {
    // 获取低阈值多边形的外壳
    const lowerShell = lowerPolygon[0];
    
    // 查找包含在这个低阈值多边形内的高阈值多边形
    const containedUpperPolygons = upperPolygons.filter(upperPolygon => {
      // 获取高阈值多边形的一个内部点
      const upperPoint = getPointInsideRing(upperPolygon[0]);
      
      // 检查这个点是否在低阈值多边形内
      return isPointInPolygon(upperPoint, lowerShell);
    });
    
    // 如果没有包含的高阈值多边形，直接添加低阈值多边形
    if (containedUpperPolygons.length === 0) {
      resultPolygons.push(lowerPolygon);
      continue;
    }
    
    // 创建一个新的多边形，使用低阈值多边形的外壳作为外壳
    const newPolygon = [lowerShell];
    
    // 添加低阈值多边形的洞
    for (let i = 1; i < lowerPolygon.length; i++) {
      newPolygon.push(lowerPolygon[i]);
    }
    
    // 添加包含的高阈值多边形的外壳作为洞
    for (const upperPolygon of containedUpperPolygons) {
      newPolygon.push(upperPolygon[0]);
      
      // 添加高阈值多边形的洞作为低阈值多边形的外壳
      for (let i = 1; i < upperPolygon.length; i++) {
        resultPolygons.push([upperPolygon[i]]);
      }
    }
    
    // 添加新的多边形
    resultPolygons.push(newPolygon);
  }
  
  // 添加不包含在任何低阈值多边形内的高阈值多边形
  for (const upperPolygon of upperPolygons) {
    const upperPoint = getPointInsideRing(upperPolygon[0]);
    
    let isContained = false;
    for (const lowerPolygon of lowerPolygons) {
      if (isPointInPolygon(upperPoint, lowerPolygon[0])) {
        isContained = true;
        break;
      }
    }
    
    if (!isContained) {
      resultPolygons.push(upperPolygon);
    }
  }
  
  return resultPolygons;
}

/**
 * 简化多边形
 * @param {Array<Array<Number>>} polygon 多边形
 * @param {Number} tolerance 容差
 * @returns {Array<Array<Number>>} 简化后的多边形
 */
export function simplifyPolygon(polygon, tolerance = 0.1) {
  if (!polygon || polygon.length < 3) {
    return polygon;
  }
  
  // 使用Douglas-Peucker算法简化多边形
  return douglasPeucker(polygon, tolerance);
}

/**
 * Douglas-Peucker算法
 * @param {Array<Array<Number>>} points 点数组
 * @param {Number} tolerance 容差
 * @returns {Array<Array<Number>>} 简化后的点数组
 */
function douglasPeucker(points, tolerance) {
  // 如果点数小于3，无法简化
  if (points.length < 3) {
    return points;
  }
  
  // 查找最大距离点
  let maxDistance = 0;
  let index = 0;
  
  const firstPoint = points[0];
  const lastPoint = points[points.length - 1];
  
  for (let i = 1; i < points.length - 1; i++) {
    const distance = perpendicularDistance(points[i], firstPoint, lastPoint);
    
    if (distance > maxDistance) {
      maxDistance = distance;
      index = i;
    }
  }
  
  // 如果最大距离大于容差，递归简化
  if (maxDistance > tolerance) {
    // 递归处理两部分
    const firstPart = douglasPeucker(points.slice(0, index + 1), tolerance);
    const lastPart = douglasPeucker(points.slice(index), tolerance);
    
    // 合并结果，去除重复点
    return [...firstPart.slice(0, -1), ...lastPart];
  } else {
    // 如果最大距离小于容差，只保留首尾点
    return [firstPoint, lastPoint];
  }
}

/**
 * 计算点到线段的垂直距离
 * @param {Array<Number>} point 点
 * @param {Array<Number>} lineStart 线段起点
 * @param {Array<Number>} lineEnd 线段终点
 * @returns {Number} 距离
 */
function perpendicularDistance(point, lineStart, lineEnd) {
  const x = point[0];
  const y = point[1];
  const x1 = lineStart[0];
  const y1 = lineStart[1];
  const x2 = lineEnd[0];
  const y2 = lineEnd[1];
  
  const dx = x2 - x1;
  const dy = y2 - y1;
  
  // 如果线段是一个点，直接计算点到点的距离
  if (dx === 0 && dy === 0) {
    const d1 = x - x1;
    const d2 = y - y1;
    return Math.sqrt(d1 * d1 + d2 * d2);
  }
  
  // 计算垂直距离
  const lineLengthSquared = dx * dx + dy * dy;
  const t = ((x - x1) * dx + (y - y1) * dy) / lineLengthSquared;
  
  if (t < 0) {
    // 点到起点的距离
    const d1 = x - x1;
    const d2 = y - y1;
    return Math.sqrt(d1 * d1 + d2 * d2);
  }
  
  if (t > 1) {
    // 点到终点的距离
    const d1 = x - x2;
    const d2 = y - y2;
    return Math.sqrt(d1 * d1 + d2 * d2);
  }
  
  // 点到线段的垂直距离
  const projX = x1 + t * dx;
  const projY = y1 + t * dy;
  const d1 = x - projX;
  const d2 = y - projY;
  
  return Math.sqrt(d1 * d1 + d2 * d2);
} 