# 增强版等值线与等值面生成器

这个库提供了一个基于Marching Squares算法的等值线和等值面生成器实现，参考了[openhome.cc](https://openhome.cc/Gossip/P5JS/MarchingSquares.html)系列文章的思路，增加了多项改进功能。

## 特点

1. **平滑的等值线**：使用了插值和Chaikin平滑算法生成更平滑的轮廓线，没有生硬的折点
2. **支持null值**：能够处理包含null、undefined或NaN的数据
3. **鞍点二义性解决**：解决了鞍点情况下连接方式的二义性问题
4. **边界处理**：改进的边界处理，确保完整的轮廓生成
5. **多边形拓扑**：正确处理等值面的内外环关系
6. **GeoJSON格式输出**：输出符合GeoJSON格式的数据，便于与可视化库集成

## 使用方法

### 等值线生成

```javascript
import { generateContours } from './enhanced-contour.js';

// 示例数据 - 二维数组，可以包含null值
const data = [
  [0.0, 0.2, 0.4, 0.6, 0.8],
  [0.2, 0.4, null, 0.8, 1.0],
  [0.4, 0.6, 0.8, 1.0, 1.2],
  [0.6, 0.8, 1.0, 1.2, 1.4],
  [0.8, 1.0, 1.2, 1.4, 1.6]
];

// 生成阈值为0.5的等值线
const contours = generateContours(data, 0.5);

console.log(contours);
// 输出格式: 
// {
//   type: "MultiLineString",
//   threshold: 0.5,
//   coordinates: [...] // 等值线坐标数组
// }
```

### 等值面生成

```javascript
import { generateContourBands } from './enhanced-contour.js';

// 使用相同的示例数据
const data = [ /* ... */ ];

// 生成阈值区间为[0.5, 1.0]的等值面
const bands = generateContourBands(data, 0.5, 1.0);

console.log(bands);
// 输出格式:
// {
//   type: "MultiPolygon",
//   lowerValue: 0.5,
//   upperValue: 1.0,
//   coordinates: [...] // 等值面坐标数组，包含外环和内环(洞)
// }
```

## 算法原理

### Marching Squares
这个实现基于Marching Squares算法，该算法将2D数据网格分成小单元格，并基于四个角点的值确定等值线如何穿过单元格。它将处理16种可能的情况：

```
0: □□    1: ■□    2: □■    3: ■■    4: □□    5: ■□    6: □■    7: ■■
  □□      □□      □□      □□      ■□      ■□      ■□      ■□

8: □□    9: ■□    10: □■   11: ■■   12: □□   13: ■□   14: □■   15: ■■
  □■      □■      □■      □■      ■■      ■■      ■■      ■■
```

其中黑色(■)表示高于阈值，白色(□)表示低于阈值的点。

### 鞍点处理

当遇到情况5和10(鞍点)时，连接方式有两种可能，我们通过计算单元格中心点的值来决定使用哪种连接方式：

```
情况5:        或        情况10:       或
■□           ■---□      □■           □---■
|     VS      |          |     VS      |
□---■         ■□         ■---□         □■
```

### 平滑算法

我们使用两种方法来生成更平滑的等值线：

1. **插值优化**：使用平滑插值函数而非简单的线性插值，在计算等值线穿过单元格边缘的位置时产生更自然的曲线
2. **Chaikin平滑**：对生成的轮廓应用Chaikin平滑算法，在每条线段上递归插入更多点，使轮廓更加圆滑

## 演示

参见`enhanced-example.js`文件，其中包含了使用各种数据集生成等值线和等值面的完整示例：

```javascript
// 运行所有演示
import { runAllDemos } from './enhanced-example.js';
runAllDemos();

// 或者单独运行特定演示
import { demoContours, demoContourBands } from './enhanced-example.js';
demoContours();
demoContourBands();
```

## 与原版d3-contour的区别

相比d3-contour，这个实现：

1. 代码结构更清晰，更接近文章中的实现思路
2. 生成的等值线更平滑，具有更好的视觉效果
3. 简化了API，使其更易于使用
4. 更好地处理边缘情况和特殊输入

## 参考资料

- [Marching squares（一）](https://openhome.cc/Gossip/P5JS/MarchingSquares.html)
- [Marching squares（二）](https://openhome.cc/Gossip/P5JS/MarchingSquares2.html)
- [Marching squares（三）](https://openhome.cc/Gossip/P5JS/MarchingSquares3.html)
- [Wikipedia: Marching squares](https://en.wikipedia.org/wiki/Marching_squares)
- [Chaikin's Algorithm for Curve Smoothing](https://www.cs.unc.edu/~dm/UNC/COMP258/LECTURES/Chaikins-Algorithm.pdf) 