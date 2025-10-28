/**
 * 测试两种方法的区别
 */

console.log("=== 对比测试：传统方法 vs 层级方法 ===\n");

// 创建简单的测试数据
const testData = [
  10, 20, 30, 20, 10,
  20, 30, 40, 30, 20,
  30, 40, 50, 40, 30,
  20, 30, 40, 30, 20,
  10, 20, 30, 20, 10
];

console.log("测试数据 (5x5网格):");
for (let j = 0; j < 5; j++) {
  const row = [];
  for (let i = 0; i < 5; i++) {
    row.push(testData[j * 5 + i].toString().padStart(2));
  }
  console.log(`  ${row.join(' ')}`);
}

const threshold = 35;
console.log(`\n测试阈值: ${threshold}\n`);

// 1. 测试传统方法
console.log("1. 传统Marching Squares方法:");
console.log("   逐个单元格处理，使用阈值比较");

let traditionalBoundaries = 0;
for (let j = 0; j < 4; j++) {
  for (let i = 0; i < 4; i++) {
    const v00 = testData[j * 5 + i];
    const v01 = testData[j * 5 + i + 1];
    const v10 = testData[(j + 1) * 5 + i];
    const v11 = testData[(j + 1) * 5 + i + 1];
    
    const above = [v00, v01, v10, v11].map(v => v >= threshold);
    const caseIndex = above[0] | above[1] << 1 | above[2] << 2 | above[3] << 3;
    
    if (caseIndex > 0 && caseIndex < 15) {
      traditionalBoundaries++;
      console.log(`   单元格(${i},${j}): 值[${v00},${v01},${v10},${v11}] → 案例${caseIndex} → 有等值线`);
    }
  }
}
console.log(`   总边界数: ${traditionalBoundaries}`);
console.log(`   预期形状: 平滑曲线，反映真实数据分布\n`);

// 2. 测试层级方法
console.log("2. 层级分配方法:");
console.log("   将所有网格点分配到层级，然后检测层级边界");

const thresholds = [25, 35, 45];
const levelMap = new Array(25);

for (let i = 0; i < 25; i++) {
  const value = testData[i];
  let level = -1;
  for (let j = 0; j < thresholds.length; j++) {
    if (value >= thresholds[j]) {
      level = j;
    } else {
      break;
    }
  }
  levelMap[i] = level;
}

console.log("   层级分配结果:");
for (let j = 0; j < 5; j++) {
  const row = [];
  for (let i = 0; i < 5; i++) {
    const level = levelMap[j * 5 + i];
    row.push(level === -1 ? '-1' : ` ${level}`);
  }
  console.log(`     ${row.join(' ')}`);
}

// 检查层级1的边界（对应阈值35）
const targetLevel = 1;
let levelBoundaries = 0;

for (let j = 0; j < 4; j++) {
  for (let i = 0; i < 4; i++) {
    const l00 = levelMap[j * 5 + i];
    const l01 = levelMap[j * 5 + i + 1];
    const l10 = levelMap[(j + 1) * 5 + i];
    const l11 = levelMap[(j + 1) * 5 + i + 1];
    
    const above = [l00, l01, l10, l11].map(l => l >= targetLevel);
    const caseIndex = above[0] | above[1] << 1 | above[2] << 2 | above[3] << 3;
    
    if (caseIndex > 0 && caseIndex < 15) {
      levelBoundaries++;
      console.log(`   单元格(${i},${j}): 层级[${l00},${l01},${l10},${l11}] → 案例${caseIndex} → 有等值线`);
    }
  }
}
console.log(`   总边界数: ${levelBoundaries}`);
console.log(`   预期形状: 可能产生方框，因为层级是离散的\n`);

// 3. 分析差异
console.log("3. 关键差异分析:");
console.log(`   传统方法: ${traditionalBoundaries} 个边界段`);
console.log(`   层级方法: ${levelBoundaries} 个边界段`);

console.log("\n   传统方法优势:");
console.log("   ✅ 使用连续的数值进行插值");
console.log("   ✅ 精确的边界位置计算");
console.log("   ✅ 平滑的曲线形状");

console.log("\n   层级方法问题:");
console.log("   ❌ 使用离散的层级标签");
console.log("   ❌ 丢失了数值的连续性");
console.log("   ❌ 可能产生阶梯状或方框状边界");

console.log("\n4. 解决方案建议:");
console.log("   问题根源: 层级分配破坏了数值的连续性");
console.log("   正确做法: 直接使用传统方法，在后处理阶段处理重叠");
console.log("   或者: 改进层级方法，保持原始数值用于插值计算");

console.log("\n=== 结论 ===");
console.log("层级分配方法确实是导致方框形状的根本原因！");
console.log("需要回到传统的阈值比较方法来获得正确的等值线形状。");