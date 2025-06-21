# Enhanced Contour 增强版等值线生成器

基于 Marching Squares 算法实现的增强版等值线和等值面生成器。参考了 openhome.cc 的 Marching Squares 系列文章，并在 d3-contour 基础上进行了改进。

## 特点

1. **完善的 null 值处理**：可以处理数据中的 null 或 NaN 值，避免生成错误的等值线
2. **平滑等值线**：使用平滑曲线插值和 Chaikin 算法生成更流畅的等值线
3. **鞍点二义性处理**：通过计算单元格中心值解决鞍点情况（case 5 和 case 10）的二义性问题
4. **改进的边界处理**：确保等值线在边界处正确闭合
5. **完善的拓扑关系**：正确处理等值面的内外环关系

## 最近更新

- **等值面渲染修复**：修复了等值面的渲染范围和闭合构造问题
  - 添加了数据边界轮廓作为最外层轮廓
  - 改进了内外环关系的处理逻辑
  - 确保环的正确闭合和方向
  - 处理边界外的区域，确保完整覆盖

## 主要API

### generateContours(data, threshold)

生成指定阈值的等值线。

参数：
- `data`: 二维数组，表示数据网格
- `threshold`: 数值，表示等值线的阈值

返回：
- GeoJSON格式的等值线对象，类型为 "MultiLineString"

### generateContourBands(data, lowerThreshold, upperThreshold)

生成指定阈值范围的等值面。

参数：
- `data`: 二维数组，表示数据网格
- `lowerThreshold`: 数值，表示等值面的下阈值
- `upperThreshold`: 数值，表示等值面的上阈值

返回：
- GeoJSON格式的等值面对象，类型为 "MultiPolygon"

## 使用示例

```javascript
import { generateContours, generateContourBands } from './enhanced-contour.js';

// 示例数据
const data = [
  [0.0, 0.2, 0.4, 0.6, 0.8],
  [0.2, 0.4, 0.6, 0.8, 1.0],
  [0.4, 0.6, null, 1.0, 1.2],
  [0.6, 0.8, 1.0, 1.2, 1.4],
  [0.8, 1.0, 1.2, 1.4, 1.6]
];

// 生成等值线
const contours = generateContours(data, 0.5);
console.log(`生成了${contours.coordinates.length}条等值线`);

// 生成等值面
const bands = generateContourBands(data, 0.5, 1.0);
console.log(`生成了${bands.coordinates.length}个等值面`);
```

## 与D3-Contour的区别

1. **null值处理**：d3-contour不能很好地处理null值，而enhanced-contour可以
2. **平滑曲线**：enhanced-contour生成的等值线更加平滑
3. **鞍点处理**：enhanced-contour通过计算中心值解决鞍点二义性
4. **边界处理**：enhanced-contour改进了边界处理，确保轮廓正确闭合
5. **等值面处理**：enhanced-contour正确处理等值面的内外环关系和边界区域

## 实现细节

### Marching Squares算法

Marching Squares是一种生成等值线的经典算法。它将数据网格分成一个个小方格（单元格），并根据单元格四个角点的值与阈值的关系，确定等值线如何穿过该单元格。

### 鞍点处理

当单元格的对角两个角点高于阈值，另外两个角点低于阈值时，会出现鞍点。这种情况下，等值线的连接方式存在二义性。enhanced-contour通过计算单元格中心点的值来解决这个问题。

### 平滑算法

enhanced-contour使用两种方法使等值线更平滑：
1. 在插值计算等值线交点时使用平滑曲线插值
2. 对生成的等值线应用Chaikin平滑算法

### 等值面生成

等值面是由两个阈值之间的区域形成的。enhanced-contour通过以下步骤生成等值面：
1. 生成下阈值的等值线作为外环
2. 生成上阈值的等值线作为内环（洞）
3. 确定内外环的包含关系
4. 处理边界区域，确保完整覆盖数据范围

## 更多示例

详见 `enhanced-example.js` 文件，其中包含了等值线、等值面以及多级等值面的生成和绘制示例。

## 参考资料

- [Marching squares（一）](https://openhome.cc/Gossip/P5JS/MarchingSquares.html)
- [Marching squares（二）](https://openhome.cc/Gossip/P5JS/MarchingSquares2.html)
- [Marching squares（三）](https://openhome.cc/Gossip/P5JS/MarchingSquares3.html)
- [Wikipedia: Marching squares](https://en.wikipedia.org/wiki/Marching_squares)
- [Chaikin's Algorithm for Curve Smoothing](https://www.cs.unc.edu/~dm/UNC/COMP258/LECTURES/Chaikins-Algorithm.pdf) 