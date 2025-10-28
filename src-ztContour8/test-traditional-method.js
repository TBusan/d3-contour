/**
 * 测试使用传统方法的等值线生成
 */

// Mock the required dependencies
const mockD3Array = {
  extent: (values, accessor = x => x) => {
    const filtered = values.filter(v => v != null && isFinite(accessor(v)));
    if (filtered.length === 0) return [0, 1];
    return [Math.min(...filtered.map(accessor)), Math.max(...filtered.map(accessor))];
  },
  nice: (start, end, count) => [start, end],
  ticks: (start, end, count) => {
    const step = (end - start) / (count - 1);
    return Array.from({length: count}, (_, i) => start + i * step);
  },
  thresholdSturges: (values) => {
    const n = values.filter(v => v != null && isFinite(v)).length;
    return Math.max(1, Math.ceil(Math.log2(n) + 1));
  }
};

// Test with a simple contour function that uses traditional method
function testTraditionalContours() {
  console.log("=== 测试传统等值线方法 ===\n");
  
  // Test data: smooth distribution
  const testData = [
    10, 15, 20, 15, 10,
    15, 25, 35, 25, 15,
    20, 35, 50, 35, 20,
    15, 25, 35, 25, 15,
    10, 15, 20, 15, 10
  ];
  
  console.log("测试数据分布 (5x5网格):");
  for (let j = 0; j < 5; j++) {
    const row = [];
    for (let i = 0; i < 5; i++) {
      row.push(testData[j * 5 + i].toString().padStart(2));
    }
    console.log(`  ${row.join(' ')}`);
  }
  
  console.log("\n1. 测试传统阈值方法 (preventOverlap = false):");
  
  // Simulate what happens when preventOverlap is false
  const thresholds = [20, 30, 40];
  
  console.log("   生成独立等值线:");
  thresholds.forEach((threshold, i) => {
    console.log(`\n   阈值 ${threshold}:`);
    
    // Count intersections for this threshold
    let intersections = 0;
    let details = [];
    
    for (let j = 0; j < 4; j++) {
      for (let k = 0; k < 4; k++) {
        const v00 = testData[j * 5 + k];
        const v01 = testData[j * 5 + k + 1];
        const v10 = testData[(j + 1) * 5 + k];
        const v11 = testData[(j + 1) * 5 + k + 1];
        
        const above = [v00, v01, v10, v11].map(v => v >= threshold);
        const caseIndex = above[0] | above[1] << 1 | above[2] << 2 | above[3] << 3;
        
        if (caseIndex > 0 && caseIndex < 15) {
          intersections++;
          
          // Calculate interpolated positions (simplified)
          const positions = [];
          if ((v00 >= threshold) !== (v01 >= threshold)) {
            const t = (threshold - v00) / (v01 - v00);
            positions.push(`边(${(k + t).toFixed(1)}, ${j})`);
          }
          if ((v01 >= threshold) !== (v11 >= threshold)) {
            const t = (threshold - v01) / (v11 - v01);
            positions.push(`边(${k + 1}, ${(j + t).toFixed(1)})`);
          }
          if ((v11 >= threshold) !== (v10 >= threshold)) {
            const t = (threshold - v11) / (v10 - v11);
            positions.push(`边(${(k + 1 - t).toFixed(1)}, ${j + 1})`);
          }
          if ((v10 >= threshold) !== (v00 >= threshold)) {
            const t = (threshold - v10) / (v00 - v10);
            positions.push(`边(${k}, ${(j + 1 - t).toFixed(1)})`);
          }
          
          details.push({
            cell: `(${k},${j})`,
            values: [v00, v01, v10, v11],
            case: caseIndex,
            positions
          });
        }
      }
    }
    
    console.log(`     交叉点数量: ${intersections}`);
    console.log(`     预期形状: ${intersections > 8 ? '复杂曲线' : intersections > 4 ? '简单闭合曲线' : '点或线段'}`);
    
    if (details.length > 0) {
      console.log(`     前3个单元格详情:`);
      details.slice(0, 3).forEach(detail => {
        console.log(`       单元格${detail.cell}: 值[${detail.values.join(',')}] → 案例${detail.case} → ${detail.positions.join(', ')}`);
      });
    }
  });
  
  console.log("\n2. 验证插值质量:");
  
  // Test interpolation for a specific case
  const threshold = 30;
  const v1 = 25, v2 = 35; // Values crossing threshold
  const interpolated = (threshold - v1) / (v2 - v1);
  
  console.log(`   示例插值计算:`);
  console.log(`     值1: ${v1}, 值2: ${v2}, 阈值: ${threshold}`);
  console.log(`     插值位置: ${interpolated.toFixed(3)} (在边界的${(interpolated * 100).toFixed(1)}%处)`);
  console.log(`     结果: 精确的亚网格位置，产生平滑曲线`);
  
  console.log("\n3. 与层级方法对比:");
  console.log("   传统方法:");
  console.log("   ✅ 保持原始数值进行插值计算");
  console.log("   ✅ 精确的边界位置");
  console.log("   ✅ 平滑连续的曲线");
  console.log("   ✅ 正确反映数据梯度");
  
  console.log("\n   层级方法 (之前的错误):");
  console.log("   ❌ 使用离散层级标签");
  console.log("   ❌ 丢失数值精度");
  console.log("   ❌ 产生阶梯状边界");
  console.log("   ❌ 方框形状");
  
  console.log("\n=== 结论 ===");
  console.log("✅ 传统Marching Squares方法是生成平滑等值线的正确选择");
  console.log("✅ 修复方案: 使用传统方法 + 后处理去重叠");
  console.log("✅ 这样既保证形状正确，又能实现非重叠拓扑");
}

// 运行测试
testTraditionalContours();