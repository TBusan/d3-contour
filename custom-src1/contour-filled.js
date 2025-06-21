/**
 * 填充模式等值线模块
 * 基于Plotly.js的实现
 */

import { marchingSquares } from './algorithms/marching-squares.js';
import { createPolygons } from './utils/polygon-builder.js';
import { createLookupTable } from './utils/color-scale.js';
import { mergePolygons } from './utils/polygon-utils.js';

/**
 * 生成填充模式等值线（类似于Plotly.js的实现）
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Array<Number>|Number} thresholds 阈值数组或单一阈值
 * @param {Object} options 配置选项
 * @returns {Array<Object>|Object} GeoJSON格式的填充等值线
 */
export function generateFilledContours(data, thresholds, options = {}) {
  if (!data || data.length === 0 || data[0].length === 0) {
    return { type: "MultiPolygon", coordinates: [], fills: [] };
  }
  
  // 默认配置
  const defaultOptions = {
    smooth: true,           // 是否平滑等值线
    smoothFactor: 0.25,     // 平滑因子
    connectEnds: true,      // 是否连接端点
    colorScale: 'viridis',  // 颜色比例尺
    fillOpacity: 0.7,       // 填充透明度
    showLines: true,        // 是否显示轮廓线
    lineWidth: 0.5,         // 轮廓线宽度
    lineColor: 'rgba(0, 0, 0, 0.3)', // 轮廓线颜色
    fillMode: 'toself'      // 填充模式: 'toself', 'tonext', 'tozeroy'
  };
  
  const config = { ...defaultOptions, ...options };
  
  // 如果传入的是单一阈值，转换为数组
  const isArray = Array.isArray(thresholds);
  const thresholdArray = isArray ? [...thresholds] : [thresholds];
  
  // 确保阈值数组已排序
  const sortedThresholds = thresholdArray.sort((a, b) => a - b);
  
  // 获取数据范围
  let minValue = Infinity;
  let maxValue = -Infinity;
  
  const width = data[0].length;
  const height = data.length;
  
  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      const value = data[y][x];
      if (value != null) {
        minValue = Math.min(minValue, value);
        maxValue = Math.max(maxValue, value);
      }
    }
  }
  
  // 创建颜色查找表
  const colorLookup = createLookupTable(config.colorScale);
  
  // 生成所有等值线
  const contourLines = sortedThresholds.map(threshold => {
    return {
      threshold,
      lines: marchingSquares(data, threshold, config)
    };
  });
  
  // 根据填充模式生成填充区域
  let fillRegions = [];
  
  if (config.fillMode === 'toself') {
    // 每个等值线自成一个填充区域
    fillRegions = contourLines.map(({ threshold, lines }) => {
      const polygons = createPolygons(lines, data, threshold, width, height);
      const normalizedValue = (threshold - minValue) / (maxValue - minValue);
      const color = typeof config.colorScale === 'function' 
        ? config.colorScale(normalizedValue)
        : colorLookup(normalizedValue);
        
      return {
        threshold,
        polygons,
        color,
        opacity: config.fillOpacity
      };
    });
  } else if (config.fillMode === 'tonext') {
    // 相邻等值线之间形成填充区域
    for (let i = 0; i < contourLines.length - 1; i++) {
      const lowerThreshold = contourLines[i].threshold;
      const upperThreshold = contourLines[i + 1].threshold;
      const lowerLines = contourLines[i].lines;
      const upperLines = contourLines[i + 1].lines;
      
      // 创建两个等值线之间的填充区域
      const polygons = createPolygonsBetweenContours(
        lowerLines, upperLines, data, lowerThreshold, upperThreshold, width, height
      );
      
      // 使用中间值来确定颜色
      const middleThreshold = (lowerThreshold + upperThreshold) / 2;
      const normalizedValue = (middleThreshold - minValue) / (maxValue - minValue);
      const color = typeof config.colorScale === 'function' 
        ? config.colorScale(normalizedValue)
        : colorLookup(normalizedValue);
      
      fillRegions.push({
        lowerThreshold,
        upperThreshold,
        polygons,
        color,
        opacity: config.fillOpacity
      });
    }
  } else if (config.fillMode === 'tozeroy') {
    // 每个等值线与y=0之间形成填充区域
    // 这需要额外的处理，类似于Plotly.js的实现
    // ...此处省略具体实现...
  }
  
  // 构建最终结果
  const result = {
    type: "MultiPolygon",
    coordinates: [],
    fills: []
  };
  
  // 添加所有填充区域
  fillRegions.forEach(region => {
    result.coordinates.push(...region.polygons);
    result.fills.push({
      threshold: region.threshold,
      lowerThreshold: region.lowerThreshold,
      upperThreshold: region.upperThreshold,
      color: region.color,
      opacity: region.opacity
    });
  });
  
  // 如果需要显示轮廓线，添加轮廓线信息
  if (config.showLines) {
    result.lines = contourLines.map(({ threshold, lines }) => {
      return {
        threshold,
        coordinates: lines,
        color: config.lineColor,
        width: config.lineWidth
      };
    });
  }
  
  return result;
}

/**
 * 创建两个等值线之间的填充多边形
 * @param {Array} lowerLines 低阈值等值线
 * @param {Array} upperLines 高阈值等值线
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Number} lowerThreshold 低阈值
 * @param {Number} upperThreshold 高阈值
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Array} 多边形数组
 */
function createPolygonsBetweenContours(lowerLines, upperLines, data, lowerThreshold, upperThreshold, width, height) {
  // 创建低阈值和高阈值的多边形
  const lowerPolygons = createPolygons(lowerLines, data, lowerThreshold, width, height);
  const upperPolygons = createPolygons(upperLines, data, upperThreshold, width, height);
  
  // 合并多边形，创建填充区域
  return mergePolygons(lowerPolygons, upperPolygons, data, lowerThreshold, upperThreshold);
} 