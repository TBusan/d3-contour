/**
 * 增强版等值线与等值面生成器示例
 */

import { generateContours, generateContourBands } from './enhanced-contour.js';

// 示例数据：包含一些null值的二维数组
const exampleData = [
  [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
  [0.2, 0.4, 0.6, 0.8, 1.0, 1.2],
  [0.4, 0.6, null, 1.0, 1.2, 1.4],
  [0.6, 0.8, 1.0, 1.2, 1.4, 1.6],
  [0.8, 1.0, 1.2, null, 1.6, 1.8],
  [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
];

/**
 * 生成一个高斯分布的测试数据集，包含一些随机的null值
 * @param {Number} width 数据宽度
 * @param {Number} height 数据高度
 * @param {Number} nullProbability null值的概率
 * @returns {Array<Array<Number>>} 二维数据数组
 */
function generateTestData(width, height, nullProbability = 0.1) {
  const result = [];
  
  for (let y = 0; y < height; y++) {
    const row = [];
    for (let x = 0; x < width; x++) {
      if (Math.random() < nullProbability) {
        row.push(null);
        continue;
      }
      
      // 生成高斯分布值
      const cx = width / 2;
      const cy = height / 2;
      const sigma = Math.min(width, height) / 3;
      
      // 计算到中心的距离
      const dx = x - cx;
      const dy = y - cy;
      const distance = Math.sqrt(dx * dx + dy * dy);
      
      // 生成高斯分布值
      const value = Math.exp(-(distance * distance) / (2 * sigma * sigma));
      
      row.push(value);
    }
    result.push(row);
  }
  
  return result;
}

/**
 * 演示等值线生成
 */
function demoContours() {
  console.log("=== 等值线生成示例 ===");
  
  // 使用小示例数据
  console.log("使用示例数据：");
  const contour1 = generateContours(exampleData, 1.0);
  console.log(`阈值 1.0 生成了 ${contour1.coordinates.length} 条等值线`);
  
  // 使用生成的数据
  console.log("\n使用生成数据：");
  const largeData = generateTestData(20, 20, 0.05);
  
  // 生成多个等值线
  const thresholds = [0.2, 0.4, 0.6, 0.8];
  
  thresholds.forEach(threshold => {
    const contour = generateContours(largeData, threshold);
    console.log(`阈值 ${threshold} 生成了 ${contour.coordinates.length} 条等值线`);
    
    // 输出第一条等值线的点数
    if (contour.coordinates.length > 0) {
      console.log(`  - 第一条等值线有 ${contour.coordinates[0].length} 个点`);
    }
  });
}

/**
 * 演示等值面生成
 */
function demoContourBands() {
  console.log("\n=== 等值面生成示例 ===");
  
  // 使用小示例数据
  console.log("使用示例数据：");
  const band1 = generateContourBands(exampleData, 0.8, 1.2);
  console.log(`阈值区间 [0.8, 1.2] 生成了 ${band1.coordinates.length} 个等值面`);
  
  // 使用生成的数据
  console.log("\n使用生成数据：");
  const largeData = generateTestData(20, 20, 0.05);
  
  // 生成多个等值面
  const thresholdPairs = [
    [0.1, 0.3],
    [0.3, 0.5],
    [0.5, 0.7],
    [0.7, 0.9]
  ];
  
  thresholdPairs.forEach(([lower, upper]) => {
    const band = generateContourBands(largeData, lower, upper);
    console.log(`阈值区间 [${lower}, ${upper}] 生成了 ${band.coordinates.length} 个等值面`);
    
    // 输出第一个等值面的信息
    if (band.coordinates.length > 0) {
      const polygon = band.coordinates[0];
      console.log(`  - 第一个等值面有 ${polygon.length} 个环（外环+内环）`);
      
      // 展示外环点数
      console.log(`  - 外环有 ${polygon[0].length} 个点`);
      
      // 展示内环数量
      if (polygon.length > 1) {
        console.log(`  - 有 ${polygon.length - 1} 个内环（洞）`);
      }
    }
  });
}

/**
 * 运行所有演示
 */
function runAllDemos() {
  demoContours();
  demoContourBands();
  
  return {
    message: "演示完成，请检查控制台输出"
  };
}

// 如果在浏览器环境中，将函数绑定到window对象
if (typeof window !== 'undefined') {
  window.runAllDemos = runAllDemos;
  window.demoContours = demoContours;
  window.demoContourBands = demoContourBands;
  window.generateTestData = generateTestData;
  console.log("演示函数已准备好，请调用 runAllDemos() 运行所有演示");
} else if (typeof require !== 'undefined') {
  // Node.js环境下直接运行
  runAllDemos();
}

export { runAllDemos, demoContours, demoContourBands, generateTestData }; 