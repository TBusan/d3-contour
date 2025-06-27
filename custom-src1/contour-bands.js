/**
 * 等值面生成模块
 * 实现填充区域
 */

import { marchingSquares } from './algorithms/marching-squares.js';
import { createPolygons } from './utils/polygon-builder.js';
import { createLookupTable } from './utils/color-scale.js';

/**
 * 生成等值面（带状区域）
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Array<Number>|Number} thresholds 阈值数组或单一阈值
 * @param {Object} options 配置选项
 * @returns {Array<Object>|Object} GeoJSON格式的等值面数组或单个等值面
 */
export function generateContourBands(data, thresholds, options = {}) {
  if (!data || data.length === 0 || data[0].length === 0) {
    return { type: "MultiPolygon", coordinates: [] };
  }
  
  // 默认配置
  const defaultOptions = {
    smooth: true,        // 是否平滑等值线
    smoothFactor: 0.25,  // 平滑因子
    connectEnds: true,   // 是否连接端点
    colorScale: null,    // 颜色比例尺
    fillOpacity: 0.7,    // 填充透明度
    showLines: true,     // 是否显示轮廓线
    clipToDataBounds: true, // 是否裁剪到数据边界
    extendToDataBounds: true // 是否扩展到数据边界
  };
  
  const config = { ...defaultOptions, ...options };
  
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
  
  // 获取数据范围
  let minValue = Infinity;
  let maxValue = -Infinity;
  
  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      const value = data[y][x];
      if (value != null) {
        minValue = Math.min(minValue, value);
        maxValue = Math.max(maxValue, value);
      }
    }
  }
  
  // 生成每个阈值的等值面
  for (let i = 0; i < sortedThresholds.length; i++) {
    const threshold = sortedThresholds[i];
    
    // 使用Marching Squares算法生成等值线
    const contourLines = marchingSquares(data, threshold, {
      ...config,
      saddleResolution: true, // 解决鞍点问题
      extendToDataBounds: config.extendToDataBounds // 确保等值线延伸到数据边界
    });
    
    // 从等值线创建多边形
    const polygons = createPolygons(contourLines, data, threshold, width, height);
    
    // 确定填充颜色
    let fillColor;
    if (config.colorScale) {
      if (typeof config.colorScale === 'function') {
        fillColor = config.colorScale(threshold);
      } else {
        // 使用阈值在数据范围内的相对位置来确定颜色
        const normalizedValue = (threshold - minValue) / (maxValue - minValue);
        fillColor = createLookupTable(config.colorScale)(normalizedValue);
      }
    } else {
      // 默认颜色
      fillColor = `rgba(70, 130, 180, ${config.fillOpacity})`;
    }
    
    results.push({
      type: "MultiPolygon",
      threshold: threshold,
      coordinates: polygons,
      fill: {
        color: fillColor,
        opacity: config.fillOpacity
      },
      stroke: config.showLines ? {
        color: 'rgba(0, 0, 0, 0.3)',
        width: 0.5
      } : null
    });
  }
  
  // 如果原始输入是单一阈值，返回单个结果，否则返回结果数组
  return isArray ? results : results[0];
} 