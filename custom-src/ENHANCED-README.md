# Enhanced Contour

增强版等值线和等值面生成库，基于Marching Squares算法。

## 特性

- 生成等值线（Contour Lines）
- 生成等值面（Contour Bands）
- 支持一次性生成多个阈值的等值线和等值面
- 处理缺失数据（null值）
- 高性能实现
- 支持复杂拓扑结构

## 安装

```bash
npm install enhanced-contour
```

## 使用方法

### 导入

```javascript
import { generateContours, generateContourBands, generateContourAndBands } from 'enhanced-contour';
```

### 生成等值线

```javascript
// 准备二维数据
const data = [
  [10, 20, 30, 40],
  [15, 25, 35, 45],
  [20, 30, 40, 50],
  [25, 35, 45, 55]
];

// 单一阈值的等值线
const singleContour = generateContours(data, 30);

// 多个阈值的等值线
const thresholds = [20, 30, 40];
const multipleContours = generateContours(data, thresholds);
```

### 生成等值面

```javascript
// 单一阈值的等值面
const singleBand = generateContourBands(data, 30);

// 多个阈值的等值面
const thresholds = [20, 30, 40];
const multipleBands = generateContourBands(data, thresholds);
```

### 同时生成等值线和等值面

```javascript
// 单一阈值
const singleResult = generateContourAndBands(data, 30);
const { contours: singleContour, bands: singleBand } = singleResult;

// 多个阈值
const thresholds = [20, 30, 40];
const multipleResult = generateContourAndBands(data, thresholds);
const { contours: multipleContours, bands: multipleBands } = multipleResult;
```

## 返回数据格式

### 单一阈值的等值线

```javascript
{
  type: "MultiLineString",
  threshold: 30,
  coordinates: [
    [[x1, y1], [x2, y2], ...], // 第一条等值线
    [[x1, y1], [x2, y2], ...], // 第二条等值线
    // ...
  ]
}
```

### 多个阈值的等值线

```javascript
[
  {
    type: "MultiLineString",
    threshold: 20,
    coordinates: [
      [[x1, y1], [x2, y2], ...], // 第一条等值线
      [[x2, y1], [x2, y2], ...], // 第二条等值线
      // ...
    ]
  },
  {
    type: "MultiLineString",
    threshold: 30,
    coordinates: [
      // ...
    ]
  },
  // ...
]
```

### 单一阈值的等值面

```javascript
{
  type: "MultiPolygon",
  threshold: 30,
  coordinates: [
    [
      [[x1, y1], [x2, y2], ...], // 外环
      [[x3, y3], [x4, y4], ...], // 内环（洞）
      // ...
    ],
    // 更多多边形
  ]
}
```

### 多个阈值的等值面

```javascript
[
  {
    type: "MultiPolygon",
    threshold: 20,
    coordinates: [
      [
        [[x1, y1], [x2, y2], ...], // 外环
        [[x3, y3], [x4, y4], ...], // 内环（洞）
        // ...
      ],
      // 更多多边形
    ]
  },
  {
    type: "MultiPolygon",
    threshold: 30,
    coordinates: [
      // ...
    ]
  },
  // ...
]
```

## 高级用法

### 生成自定义阈值的等值线

```javascript
// 计算数据范围
const values = data.flat().filter(v => v != null);
const min = Math.min(...values);
const max = Math.max(...values);

// 生成10个均匀分布的阈值
const thresholdCount = 10;
const thresholds = Array.from({ length: thresholdCount }, (_, i) => 
  min + (max - min) * i / (thresholdCount - 1)
);

// 生成等值线和等值面
const contours = generateContours(data, thresholds);
const bands = generateContourBands(data, thresholds);
```

### 渲染等值线和等值面

```javascript
// 渲染等值线
contours.forEach((contour, i) => {
  // 使用不同颜色渲染不同阈值的等值线
  const color = getColorForThreshold(contour.threshold);
  
  contour.coordinates.forEach(line => {
    // 使用Canvas或SVG渲染线条
    drawLine(line, color);
  });
});

// 渲染等值面
bands.forEach((band, i) => {
  // 使用不同颜色渲染不同阈值的等值面
  const color = getColorForThreshold(band.threshold);
  
  band.coordinates.forEach(polygon => {
    // 使用Canvas或SVG渲染多边形
    drawPolygon(polygon, color);
  });
});
```

## 性能优化

- 对于大型数据集，可以先进行降采样
- 可以使用Web Worker在后台线程中生成等值线和等值面
- 对于交互式应用，可以根据视口范围只处理可见区域的数据

## 许可证

MIT 