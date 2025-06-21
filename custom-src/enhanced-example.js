/**
 * Enhanced Contour 示例代码
 * 展示如何使用增强版等值线和等值面生成器
 */
import { generateContours, generateContourBands, generateContourAndBands } from './enhanced-contour.js';

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
 * @param {HTMLCanvasElement} canvas 画布元素
 */
function contourExample(canvas) {
  const ctx = canvas.getContext('2d');
  const width = canvas.width;
  const height = canvas.height;
  
  // 生成测试数据
  const data = generateTestData();
  
  // 计算数据范围
  let min = Infinity;
  let max = -Infinity;
  
  for (let y = 0; y < data.length; y++) {
    for (let x = 0; x < data[0].length; x++) {
      const value = data[y][x];
      if (value != null) {
        min = Math.min(min, value);
        max = Math.max(max, value);
      }
    }
  }
  
  // 生成多个阈值
  const thresholds = [];
  const thresholdCount = 10;
  for (let i = 0; i < thresholdCount; i++) {
    thresholds.push(min + (max - min) * i / (thresholdCount - 1));
  }
  
  // 一次性生成所有阈值的等值线
  const contours = generateContours(data, thresholds);
  
  // 清空画布
  ctx.clearRect(0, 0, width, height);
  
  // 绘制等值线
  ctx.lineWidth = 1;
  
  contours.forEach((contour, i) => {
    // 根据阈值设置不同的颜色
    const hue = 240 * (1 - i / (thresholdCount - 1));
    ctx.strokeStyle = `hsl(${hue}, 100%, 50%)`;
    
    contour.coordinates.forEach(line => {
      ctx.beginPath();
      
      line.forEach((point, j) => {
        // 将数据坐标映射到画布坐标
        const canvasX = point[0] / data[0].length * width;
        const canvasY = point[1] / data.length * height;
        
        if (j === 0) {
          ctx.moveTo(canvasX, canvasY);
        } else {
          ctx.lineTo(canvasX, canvasY);
        }
      });
      
      ctx.stroke();
    });
  });
}

/**
 * 等值面示例
 * @param {HTMLCanvasElement} canvas 画布元素
 */
function contourBandExample(canvas) {
  const ctx = canvas.getContext('2d');
  const width = canvas.width;
  const height = canvas.height;
  
  // 生成测试数据
  const data = generateTestData();
  
  // 计算数据范围
  let min = Infinity;
  let max = -Infinity;
  
  for (let y = 0; y < data.length; y++) {
    for (let x = 0; x < data[0].length; x++) {
      const value = data[y][x];
      if (value != null) {
        min = Math.min(min, value);
        max = Math.max(max, value);
      }
    }
  }
  
  // 生成多个阈值
  const thresholds = [];
  const thresholdCount = 10;
  for (let i = 0; i < thresholdCount; i++) {
    thresholds.push(min + (max - min) * i / (thresholdCount - 1));
  }
  
  // 一次性生成所有阈值的等值面
  const bands = generateContourBands(data, thresholds);
  
  // 清空画布
  ctx.clearRect(0, 0, width, height);
  
  // 渲染等值面
  bands.forEach((band, i) => {
    // 根据阈值设置不同的颜色
    const hue = 240 * (1 - i / (thresholdCount - 1));
    ctx.fillStyle = `hsla(${hue}, 100%, 50%, 0.5)`;
    ctx.strokeStyle = `hsla(${hue}, 100%, 30%, 0.8)`;
    
    band.coordinates.forEach(polygon => {
      // 绘制多边形
      ctx.beginPath();
      
      // 绘制外环
      const outerRing = polygon[0];
      outerRing.forEach((point, j) => {
        // 将数据坐标映射到画布坐标
        const canvasX = point[0] / data[0].length * width;
        const canvasY = point[1] / data.length * height;
        
        if (j === 0) {
          ctx.moveTo(canvasX, canvasY);
        } else {
          ctx.lineTo(canvasX, canvasY);
        }
      });
      
      // 绘制内环（洞）
      for (let r = 1; r < polygon.length; r++) {
        const innerRing = polygon[r];
        
        // 移动到内环的第一个点
        const firstPoint = innerRing[0];
        const firstX = firstPoint[0] / data[0].length * width;
        const firstY = firstPoint[1] / data.length * height;
        ctx.moveTo(firstX, firstY);
        
        // 绘制内环的其余部分
        for (let j = 1; j < innerRing.length; j++) {
          const point = innerRing[j];
          const canvasX = point[0] / data[0].length * width;
          const canvasY = point[1] / data.length * height;
          ctx.lineTo(canvasX, canvasY);
        }
      }
      
      ctx.fill();
      ctx.stroke();
    });
  });
}

/**
 * 组合等值线和等值面示例
 * @param {HTMLCanvasElement} canvas 画布元素
 */
function combinedExample(canvas) {
  const ctx = canvas.getContext('2d');
  const width = canvas.width;
  const height = canvas.height;
  
  // 生成测试数据
  const data = generateTestData();
  
  // 计算数据范围
  let min = Infinity;
  let max = -Infinity;
  
  for (let y = 0; y < data.length; y++) {
    for (let x = 0; x < data[0].length; x++) {
      const value = data[y][x];
      if (value != null) {
        min = Math.min(min, value);
        max = Math.max(max, value);
      }
    }
  }
  
  // 生成多个阈值
  const thresholds = [];
  const thresholdCount = 10;
  for (let i = 0; i < thresholdCount; i++) {
    thresholds.push(min + (max - min) * i / (thresholdCount - 1));
  }
  
  // 一次性生成所有阈值的等值线和等值面
  const { contours, bands } = generateContourAndBands(data, thresholds);
  
  // 清空画布
  ctx.clearRect(0, 0, width, height);
  
  // 先绘制等值面
  bands.forEach((band, i) => {
    // 根据阈值设置不同的颜色
    const hue = 240 * (1 - i / (thresholdCount - 1));
    ctx.fillStyle = `hsla(${hue}, 100%, 50%, 0.3)`;
    
    band.coordinates.forEach(polygon => {
      // 绘制多边形
      ctx.beginPath();
      
      // 绘制外环
      const outerRing = polygon[0];
      outerRing.forEach((point, j) => {
        // 将数据坐标映射到画布坐标
        const canvasX = point[0] / data[0].length * width;
        const canvasY = point[1] / data.length * height;
        
        if (j === 0) {
          ctx.moveTo(canvasX, canvasY);
        } else {
          ctx.lineTo(canvasX, canvasY);
        }
      });
      
      // 绘制内环（洞）
      for (let r = 1; r < polygon.length; r++) {
        const innerRing = polygon[r];
        
        // 移动到内环的第一个点
        const firstPoint = innerRing[0];
        const firstX = firstPoint[0] / data[0].length * width;
        const firstY = firstPoint[1] / data.length * height;
        ctx.moveTo(firstX, firstY);
        
        // 绘制内环的其余部分
        for (let j = 1; j < innerRing.length; j++) {
          const point = innerRing[j];
          const canvasX = point[0] / data[0].length * width;
          const canvasY = point[1] / data.length * height;
          ctx.lineTo(canvasX, canvasY);
        }
      }
      
      ctx.fill();
    });
  });
  
  // 再绘制等值线
  ctx.lineWidth = 1;
  
  contours.forEach((contour, i) => {
    // 根据阈值设置不同的颜色
    const hue = 240 * (1 - i / (thresholdCount - 1));
    ctx.strokeStyle = `hsl(${hue}, 100%, 30%)`;
    
    contour.coordinates.forEach(line => {
      ctx.beginPath();
      
      line.forEach((point, j) => {
        // 将数据坐标映射到画布坐标
        const canvasX = point[0] / data[0].length * width;
        const canvasY = point[1] / data.length * height;
        
        if (j === 0) {
          ctx.moveTo(canvasX, canvasY);
        } else {
          ctx.lineTo(canvasX, canvasY);
        }
      });
      
      ctx.stroke();
    });
  });
}

// 导出示例函数
export { contourExample, contourBandExample, combinedExample, generateTestData }; 