# Enhanced Contour

一个基于Marching Squares算法实现的等值线和等值面渲染库，支持填充模式等值线渲染。

## 功能特点

- **等值线生成**：基于Marching Squares算法生成等值线
- **等值面生成**：生成填充的等值面（带状区域）
- **填充模式等值线**：类似于Plotly.js的填充模式等值线
- **多种填充模式**：支持'toself'、'tonext'和'tozeroy'等填充模式
- **平滑处理**：支持等值线和等值面的平滑处理
- **颜色映射**：内置多种颜色方案，支持自定义颜色映射

## 安装

```bash
# 使用npm安装
npm install enhanced-contour

# 或使用yarn
yarn add enhanced-contour
```

## 基本用法

### 生成等值线

```javascript
import { generateContours } from 'enhanced-contour';

// 创建二维数据数组
const data = [
  [10, 20, 30, 40],
  [15, 25, 35, 45],
  [20, 30, 40, 50],
  [25, 35, 45, 55]
];

// 定义阈值
const thresholds = [20, 30, 40];

// 生成等值线
const contours = generateContours(data, thresholds, {
  smooth: true,
  smoothFactor: 0.25
});

// 结果是GeoJSON格式的等值线数组
console.log(contours);
```

### 生成等值面

```javascript
import { generateContourBands } from 'enhanced-contour';

// 生成等值面
const bands = generateContourBands(data, thresholds, {
  colorScale: 'viridis',
  fillOpacity: 0.7,
  showLines: true
});

// 结果是GeoJSON格式的等值面数组
console.log(bands);
```

### 生成填充模式等值线

```javascript
import { generateFilledContours } from 'enhanced-contour';

// 生成填充模式等值线
const filled = generateFilledContours(data, thresholds, {
  colorScale: 'plasma',
  fillMode: 'tonext',
  showLines: true
});

// 结果包含填充区域和轮廓线
console.log(filled);
```

## API参考

### generateContours(data, thresholds, options)

生成等值线。

**参数**：
- `data`：二维数据数组
- `thresholds`：阈值数组或单一阈值
- `options`：配置选项
  - `smooth`：是否平滑等值线（默认：true）
  - `smoothFactor`：平滑因子（默认：0.25）
  - `connectEnds`：是否连接端点（默认：true）
  - `fill`：是否填充（默认：false）
  - `colorScale`：颜色比例尺（默认：null）

**返回值**：GeoJSON格式的等值线数组或单个等值线。

### generateContourBands(data, thresholds, options)

生成等值面（带状区域）。

**参数**：
- `data`：二维数据数组
- `thresholds`：阈值数组或单一阈值
- `options`：配置选项
  - `smooth`：是否平滑等值线（默认：true）
  - `smoothFactor`：平滑因子（默认：0.25）
  - `connectEnds`：是否连接端点（默认：true）
  - `colorScale`：颜色比例尺（默认：null）
  - `fillOpacity`：填充透明度（默认：0.7）
  - `showLines`：是否显示轮廓线（默认：true）

**返回值**：GeoJSON格式的等值面数组或单个等值面。

### generateFilledContours(data, thresholds, options)

生成填充模式等值线。

**参数**：
- `data`：二维数据数组
- `thresholds`：阈值数组或单一阈值
- `options`：配置选项
  - `smooth`：是否平滑等值线（默认：true）
  - `smoothFactor`：平滑因子（默认：0.25）
  - `connectEnds`：是否连接端点（默认：true）
  - `colorScale`：颜色比例尺（默认：'viridis'）
  - `fillOpacity`：填充透明度（默认：0.7）
  - `showLines`：是否显示轮廓线（默认：true）
  - `lineWidth`：轮廓线宽度（默认：0.5）
  - `lineColor`：轮廓线颜色（默认：'rgba(0, 0, 0, 0.3)'）
  - `fillMode`：填充模式（默认：'toself'）
    - `'toself'`：每个等值线自成一个填充区域
    - `'tonext'`：相邻等值线之间形成填充区域
    - `'tozeroy'`：每个等值线与y=0之间形成填充区域

**返回值**：包含填充区域和轮廓线的对象。

## 颜色方案

内置的颜色方案包括：

- `'viridis'`：Matplotlib的Viridis方案
- `'plasma'`：Matplotlib的Plasma方案
- `'inferno'`：Matplotlib的Inferno方案
- `'magma'`：Matplotlib的Magma方案
- `'rainbow'`：彩虹色方案
- `'blueRed'`：蓝色到红色
- `'hot'`：热力图方案
- `'cool'`：冷色调方案
- `'bgr'`：蓝绿红方案

## 示例

查看`examples`目录中的示例代码，了解更多用法。

## 浏览器兼容性

该库使用ES6+语法，建议在现代浏览器中使用。如果需要支持旧版浏览器，请使用Babel等工具进行转译。

## 许可证

MIT 