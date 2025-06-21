/**
 * 等值线生成模块
 * 基于Marching Squares算法实现
 */

import { marchingSquares } from './algorithms/marching-squares.js';
import { smoothContours } from './utils/smoothing.js';

/**
 * 生成等值线
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Array<Number>|Number} thresholds 阈值数组或单一阈值
 * @param {Object} options 配置选项
 * @returns {Array<Object>|Object} GeoJSON格式的等值线数组或单个等值线
 */
export function generateContours(data, thresholds, options = {}) {
  if (!data || data.length === 0 || data[0].length === 0) {
    return { type: "MultiLineString", coordinates: [] };
  }
  
  // 默认配置
  const defaultOptions = {
    smooth: true,        // 是否平滑等值线
    smoothFactor: 0.25,  // 平滑因子
    connectEnds: true,   // 是否连接端点
    fill: false,         // 是否填充
    colorScale: null     // 颜色比例尺
  };
  
  const config = { ...defaultOptions, ...options };
  
  // 如果传入的是单一阈值，转换为数组
  const isArray = Array.isArray(thresholds);
  const thresholdArray = isArray ? thresholds : [thresholds];
  
  // 为每个阈值生成等值线
  const results = thresholdArray.map(threshold => {
    // 使用Marching Squares算法生成等值线
    const contourLines = marchingSquares(data, threshold, config);
    
    // 平滑等值线（如果需要）
    const finalContours = config.smooth 
      ? smoothContours(contourLines, config.smoothFactor)
      : contourLines;
    
    // 创建结果对象
    const result = {
      type: "MultiLineString",
      threshold: threshold,
      coordinates: finalContours
    };
    
    // 如果需要填充，添加填充信息
    if (config.fill && config.colorScale) {
      const color = typeof config.colorScale === 'function' 
        ? config.colorScale(threshold)
        : config.colorScale;
      
      result.fill = {
        color: color,
        opacity: options.fillOpacity || 0.5
      };
    }
    
    return result;
  });
  
  // 如果原始输入是单一阈值，返回单个结果，否则返回结果数组
  return isArray ? results : results[0];
} 