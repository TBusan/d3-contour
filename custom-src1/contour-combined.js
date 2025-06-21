/**
 * 组合模块，用于同时生成等值线和等值面
 */

import { generateContours } from './contour-lines.js';
import { generateContourBands } from './contour-bands.js';
import { generateFilledContours } from './contour-filled.js';

/**
 * 同时生成等值线和等值面
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Array<Number>|Number} thresholds 阈值数组或单一阈值
 * @param {Object} options 配置选项
 * @returns {Object} 包含等值线和等值面的对象
 */
export function generateContourAndBands(data, thresholds, options = {}) {
  const contours = generateContours(data, thresholds, options);
  const bands = generateContourBands(data, thresholds, options);
  
  return {
    contours,
    bands
  };
}

/**
 * 生成完整的等值线可视化（包括填充）
 * @param {Array<Array<Number>>} data 二维数据数组
 * @param {Array<Number>|Number} thresholds 阈值数组或单一阈值
 * @param {Object} options 配置选项
 * @returns {Object} 完整的等值线可视化对象
 */
export function generateCompleteContours(data, thresholds, options = {}) {
  // 默认配置
  const defaultOptions = {
    renderMode: 'filled', // 'lines', 'bands', 'filled', 'combined'
    colorScale: 'viridis',
    fillOpacity: 0.7,
    showLines: true,
    lineWidth: 0.5,
    lineColor: 'rgba(0, 0, 0, 0.3)',
    fillMode: 'tonext'
  };
  
  const config = { ...defaultOptions, ...options };
  
  // 根据渲染模式选择不同的生成函数
  switch (config.renderMode) {
    case 'lines':
      return { contours: generateContours(data, thresholds, config) };
    
    case 'bands':
      return { bands: generateContourBands(data, thresholds, config) };
    
    case 'filled':
      return { filled: generateFilledContours(data, thresholds, config) };
    
    case 'combined':
    default:
      var contours = generateContours(data, thresholds, config);
      var filled = generateFilledContours(data, thresholds, {
        ...config,
        showLines: false // 避免重复的线条
      });
      
      return {
        contours,
        filled
      };
  }
} 