/**
 * Enhanced Contour 示例代码
 * 展示如何使用增强版等值线和等值面生成器
 */
import { generateContours, generateContourBands } from './enhanced-contour.js';

/**
 * 生成测试数据
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @returns {Array<Array<Number>>} 二维数据数组
 */
function generateTestData(width = 100, height = 100) {
  const data = [];
  for (let y = 0; y < height; y++) {
    const row = [];
    for (let x = 0; x < width; x++) {
      // 生成几个高斯分布的山峰
      const value = 
        100 * Math.exp(-((x - 30) ** 2 + (y - 30) ** 2) / 400) +
        50 * Math.exp(-((x - 70) ** 2 + (y - 70) ** 2) / 300) +
        80 * Math.exp(-((x - 50) ** 2 + (y - 50) ** 2) / 600);
      
      // 随机添加一些null值，测试null值处理
      row.push(Math.random() < 0.05 ? null : value);
    }
    data.push(row);
  }
  return data;
}

/**
 * 等值线示例
 * @returns {Object} 等值线结果
 */
function contourExample() {
  // 生成测试数据
  const data = generateTestData();
  
  // 设置等值线阈值
  const threshold = 50;
  
  // 生成等值线
  const contours = generateContours(data, threshold);
  
  console.log(`生成了${contours.coordinates.length}条等值线，阈值为${threshold}`);
  
  return contours;
}

/**
 * 等值面示例
 * @returns {Object} 等值面结果
 */
function contourBandExample() {
  // 生成测试数据
  const data = generateTestData();
  
  // 设置等值面的上下阈值
  const lowerThreshold = 30;
  const upperThreshold = 70;
  
  // 生成等值面
  const bands = generateContourBands(data, lowerThreshold, upperThreshold);
  
  console.log(`生成了${bands.coordinates.length}个等值面，阈值范围为[${lowerThreshold}, ${upperThreshold}]`);
  
  return bands;
}

/**
 * 多级等值面示例
 * @returns {Array<Object>} 多个等值面结果
 */
function multiLevelContourBandsExample() {
  // 生成测试数据
  const data = generateTestData();
  
  // 设置多级阈值
  const thresholds = [0, 20, 40, 60, 80, 100];
  
  // 生成多级等值面
  const bands = [];
  for (let i = 0; i < thresholds.length - 1; i++) {
    const lowerThreshold = thresholds[i];
    const upperThreshold = thresholds[i + 1];
    
    const band = generateContourBands(data, lowerThreshold, upperThreshold);
    bands.push(band);
    
    console.log(`生成了${band.coordinates.length}个等值面，阈值范围为[${lowerThreshold}, ${upperThreshold}]`);
  }
  
  return bands;
}

/**
 * 在Canvas上绘制等值线
 * @param {CanvasRenderingContext2D} ctx Canvas上下文
 * @param {Object} contours 等值线数据
 * @param {Number} scale 缩放比例
 * @param {String} color 线条颜色
 */
function drawContours(ctx, contours, scale = 1, color = 'blue') {
  ctx.strokeStyle = color;
  ctx.lineWidth = 1.5;
  
  contours.coordinates.forEach(contour => {
    ctx.beginPath();
    
    contour.forEach((point, i) => {
      const [x, y] = point;
      if (i === 0) {
        ctx.moveTo(x * scale, y * scale);
      } else {
        ctx.lineTo(x * scale, y * scale);
      }
    });
    
    ctx.stroke();
  });
}

/**
 * 在Canvas上绘制等值面
 * @param {CanvasRenderingContext2D} ctx Canvas上下文
 * @param {Object} bands 等值面数据
 * @param {Number} scale 缩放比例
 * @param {String} fillColor 填充颜色
 * @param {String} strokeColor 边框颜色
 */
function drawContourBands(ctx, bands, scale = 1, fillColor = 'rgba(0, 0, 255, 0.2)', strokeColor = 'rgba(0, 0, 255, 0.5)') {
  ctx.fillStyle = fillColor;
  ctx.strokeStyle = strokeColor;
  ctx.lineWidth = 1;
  
  bands.coordinates.forEach(polygon => {
    ctx.beginPath();
    
    // 绘制外环
    const outerRing = polygon[0];
    outerRing.forEach((point, i) => {
      const [x, y] = point;
      if (i === 0) {
        ctx.moveTo(x * scale, y * scale);
      } else {
        ctx.lineTo(x * scale, y * scale);
      }
    });
    
    // 绘制内环（洞）
    for (let i = 1; i < polygon.length; i++) {
      const hole = polygon[i];
      
      // 移动到洞的起点
      if (hole.length > 0) {
        ctx.moveTo(hole[0][0] * scale, hole[0][1] * scale);
        
        // 绘制洞的轮廓
        for (let j = 1; j < hole.length; j++) {
          ctx.lineTo(hole[j][0] * scale, hole[j][1] * scale);
        }
      }
    }
    
    ctx.fill();
    ctx.stroke();
  });
}

// 导出示例函数
export {
  generateTestData,
  contourExample,
  contourBandExample,
  multiLevelContourBandsExample,
  drawContours,
  drawContourBands
}; 