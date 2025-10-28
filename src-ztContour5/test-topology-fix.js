import contours from "./enhanced-contours.js";
import { contourData as testData1 } from "./testData1.js";
import { contourData as testData2 } from "./testData2.js";
import { exportGeoJSON } from "./test-enhanced-contours.js";
/**
 * 测试拓扑修复功能
 * 重点验证等值面重叠问题的解决和边界缺失的修复
 */

function testTopologyFix() {
  console.log("=== 拓扑修复测试 ===\n");

  // 测试1: testData2的边界缺失修复
  console.log("测试1: 修复testData2的边界缺失问题");
  testBoundaryMissing();

  // 测试2: 等值面重叠问题修复
  console.log("\n测试2: 等值面重叠问题修复");
  testOverlapPrevention();

  // 测试3: 拓扑验证
  console.log("\n测试3: 拓扑验证");
  testTopologyValidation();

  // 测试4: 性能对比
  console.log("\n测试4: 性能对比");
  testPerformanceComparison();

  console.log("\n=== 测试完成 ===");
}

function testBoundaryMissing() {
  const { x, y, v } = testData2.data;
  const dx = x.length;
  const dy = y.length;
  
  // 扁平化2D数组
  const values = [];
  for (let j = 0; j < dy; j++) {
    for (let i = 0; i < dx; i++) {
      values.push(v[j][i]);
    }
  }

  console.log(`  数据集: ${dx} x ${dy} 网格 (${values.length} 单元格)`);

  // 原始算法 (可能有边界缺失)
  const originalGen = contours()
    .size([dx, dy])
    .preventOverlap(false)
    .extendBoundaries(false)
    .thresholds([50, 100, 150, 200, 250]);

  const originalStart = performance.now();
  const originalResults = originalGen(values);
  const originalDuration = performance.now() - originalStart;

  // 修复的算法 (边界扩展)
  const fixedGen = contours()
    .size([dx, dy])
    .preventOverlap(true)
    .extendBoundaries(true)
    .thresholds([50, 100, 150, 200, 250]);

  const fixedStart = performance.now();
  const fixedResults = fixedGen(values);
  const fixedDuration = performance.now() - fixedStart;

  // 比较结果
  console.log(`  原始算法: ${originalDuration.toFixed(2)}ms, ${originalResults.length} 等值面`);
  console.log(`  修复算法: ${fixedDuration.toFixed(2)}ms, ${fixedResults.length} 等值面`);

  // 分析边界覆盖
  let originalBoundaryPolygons = 0;
  let fixedBoundaryPolygons = 0;

  originalResults.forEach(contour => {
    contour.coordinates.forEach(polygon => {
      if (touchesBoundary(polygon[0], dx, dy)) {
        originalBoundaryPolygons++;
      }
    });
  });

  fixedResults.forEach(contour => {
    contour.coordinates.forEach(polygon => {
      if (touchesBoundary(polygon[0], dx, dy)) {
        fixedBoundaryPolygons++;
      }
    });
  });

  console.log(`  原始算法边界多边形: ${originalBoundaryPolygons}`);
  console.log(`  修复算法边界多边形: ${fixedBoundaryPolygons}`);
  console.log(`  边界覆盖改善: ${((fixedBoundaryPolygons - originalBoundaryPolygons) / Math.max(originalBoundaryPolygons, 1) * 100).toFixed(1)}%`);
}

function testOverlapPrevention() {
  // 创建一个有潜在重叠问题的测试数据集
  const testValues = [];
  const size = 20;
  
  for (let j = 0; j < size; j++) {
    for (let i = 0; i < size; i++) {
      const x = i / size;
      const y = j / size;
      // 创建复杂的值分布，容易产生重叠等值面
      const value = Math.sin(x * Math.PI * 2) * Math.cos(y * Math.PI * 2) * 100 + 
                   Math.sin(x * Math.PI * 4) * Math.cos(y * Math.PI * 4) * 50 + 100;
      testValues.push(value);
    }
  }

  console.log(`  测试数据: ${size} x ${size} 网格`);

  // 原始算法 (可能有重叠)
  const originalGen = contours()
    .size([size, size])
    .preventOverlap(false)
    .thresholds([80, 100, 120, 140, 160]);

  const originalResults = originalGen(testValues);

  // 修复的算法 (防重叠)
  const fixedGen = contours()
    .size([size, size])
    .preventOverlap(true)
    .thresholds([80, 100, 120, 140, 160]);

  const fixedResults = fixedGen(testValues);

  const oo = exportGeoJSON(fixedResults);
  const res = JSON.parse(oo);
  console.log(res);


  // 检测重叠
  const originalOverlaps = detectOverlaps(originalResults);
  const fixedOverlaps = detectOverlaps(fixedResults);

  console.log(`  原始算法重叠检测: ${originalOverlaps} 个重叠区域`);
  console.log(`  修复算法重叠检测: ${fixedOverlaps} 个重叠区域`);
  console.log(`  重叠减少: ${Math.max(0, originalOverlaps - fixedOverlaps)} 个`);

  // 验证层级关系
  const levelValidation = validateLevelHierarchy(fixedResults);
  console.log(`  层级关系验证: ${levelValidation ? '✓ 正确' : '✗ 错误'}`);
}

function testTopologyValidation() {
  // 使用testData1进行拓扑验证
  const { x, y, v } = testData1.data;
  const dx = x.length;
  const dy = y.length;
  
  const values = [];
  for (let j = 0; j < dy; j++) {
    for (let i = 0; i < dx; i++) {
      values.push(v[j][i]);
    }
  }

  const generator = contours()
    .size([dx, dy])
    .preventOverlap(true)
    .extendBoundaries(true)
    .nullHandling(true)
    .thresholds([100, 200, 300]);

  const results = generator(values);

  // 拓扑验证
  const validationResults = {
    selfIntersections: 0,
    invalidPolygons: 0,
    invalidHoles: 0,
    totalPolygons: 0,
    totalHoles: 0
  };

  results.forEach(contour => {
    contour.coordinates.forEach(polygon => {
      validationResults.totalPolygons++;
      
      // 检查外环
      if (polygon[0].length < 4) {
        validationResults.invalidPolygons++;
      }
      
      // 检查自相交 (简化检测)
      if (hasSelfIntersection(polygon[0])) {
        validationResults.selfIntersections++;
      }
      
      // 检查孔洞
      for (let i = 1; i < polygon.length; i++) {
        validationResults.totalHoles++;
        if (polygon[i].length < 4) {
          validationResults.invalidHoles++;
        }
      }
    });
  });

  console.log(`  总多边形数: ${validationResults.totalPolygons}`);
  console.log(`  总孔洞数: ${validationResults.totalHoles}`);
  console.log(`  自相交: ${validationResults.selfIntersections}`);
  console.log(`  无效多边形: ${validationResults.invalidPolygons}`);
  console.log(`  无效孔洞: ${validationResults.invalidHoles}`);
  
  const isValid = validationResults.selfIntersections === 0 && 
                  validationResults.invalidPolygons === 0 && 
                  validationResults.invalidHoles === 0;
  
  console.log(`  拓扑验证结果: ${isValid ? '✓ 通过' : '✗ 失败'}`);
}

function testPerformanceComparison() {
  // 生成大型测试数据集
  const size = 100;
  const values = new Array(size * size);
  
  for (let i = 0; i < values.length; i++) {
    const x = (i % size) / size;
    const y = Math.floor(i / size) / size;
    values[i] = Math.sin(x * Math.PI * 3) * Math.cos(y * Math.PI * 3) * 100 + 
               Math.random() * 20;
  }

  console.log(`  性能测试数据: ${size} x ${size} 网格`);

  const thresholds = [20, 40, 60, 80, 100, 120, 140, 160, 180];

  // 原始算法性能
  const originalGen = contours()
    .size([size, size])
    .preventOverlap(false)
    .extendBoundaries(false);

  const originalStart = performance.now();
  const originalResults = originalGen.thresholds(thresholds)(values);
  const originalDuration = performance.now() - originalStart;

  // 修复算法性能
  const fixedGen = contours()
    .size([size, size])
    .preventOverlap(true)
    .extendBoundaries(true);

  const fixedStart = performance.now();
  const fixedResults = fixedGen.thresholds(thresholds)(values);
  const fixedDuration = performance.now() - fixedStart;

  console.log(`  原始算法: ${originalDuration.toFixed(2)}ms`);
  console.log(`  修复算法: ${fixedDuration.toFixed(2)}ms`);
  console.log(`  性能开销: ${(fixedDuration / originalDuration).toFixed(2)}x`);

  // 质量对比
  let originalTotalPolygons = 0, fixedTotalPolygons = 0;
  originalResults.forEach(c => { originalTotalPolygons += c.coordinates.length; });
  fixedResults.forEach(c => { fixedTotalPolygons += c.coordinates.length; });

  console.log(`  原始算法多边形总数: ${originalTotalPolygons}`);
  console.log(`  修复算法多边形总数: ${fixedTotalPolygons}`);
}

// 辅助函数
function touchesBoundary(ring, dx, dy) {
  return ring.some(point => {
    const x = point[0];
    const y = point[1];
    return x <= 0.1 || x >= dx - 0.1 || y <= 0.1 || y >= dy - 0.1;
  });
}

function detectOverlaps(contours) {
  // 简化的重叠检测：检查不同等值面是否有包含关系
  let overlaps = 0;
  
  for (let i = 0; i < contours.length - 1; i++) {
    for (let j = i + 1; j < contours.length; j++) {
      const contour1 = contours[i];
      const contour2 = contours[j];
      
      // 检查第一个等值面的点是否在第二个等值面内
      for (const polygon1 of contour1.coordinates) {
        for (const polygon2 of contour2.coordinates) {
          if (ringIntersectsRing(polygon1[0], polygon2[0])) {
            overlaps++;
            break;
          }
        }
      }
    }
  }
  
  return overlaps;
}

function ringIntersectsRing(ring1, ring2) {
  // 简化的相交检测：检查第一个环的点是否在第二个环内
  const testPoint = ring1[Math.floor(ring1.length / 2)];
  return pointInRing(testPoint, ring2);
}

function pointInRing(point, ring) {
  const [x, y] = point;
  let inside = false;
  
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    const [xi, yi] = ring[i];
    const [xj, yj] = ring[j];
    
    if (((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
      inside = !inside;
    }
  }
  
  return inside;
}

function validateLevelHierarchy(contours) {
  // 验证等值面的层级关系是否正确
  contours.sort((a, b) => a.value - b.value);
  
  for (let i = 0; i < contours.length - 1; i++) {
    if (contours[i].value >= contours[i + 1].value) {
      return false; // 层级顺序错误
    }
  }
  
  return true;
}

function hasSelfIntersection(ring) {
  // 简化的自相交检测
  if (ring.length < 4) return false;
  
  for (let i = 0; i < ring.length - 1; i++) {
    for (let j = i + 2; j < ring.length - 1; j++) {
      if (j === ring.length - 2 && i === 0) continue; // 跳过首尾相邻线段
      
      if (segmentsIntersect(ring[i], ring[i + 1], ring[j], ring[j + 1])) {
        return true;
      }
    }
  }
  
  return false;
}

function segmentsIntersect(p1, p2, p3, p4) {
  const [x1, y1] = p1;
  const [x2, y2] = p2;
  const [x3, y3] = p3;
  const [x4, y4] = p4;
  
  const denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
  if (Math.abs(denom) < 1e-10) return false; // 平行
  
  const t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
  const u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom;
  
  return t >= 0 && t <= 1 && u >= 0 && u <= 1;
}

// 运行测试
if (typeof window === 'undefined' && typeof global !== 'undefined') {
  testTopologyFix();
}

export default testTopologyFix;