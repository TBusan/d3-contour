# Enhanced Contour 等值线实现分析

## 概述

`enhanced-contour.js` 是基于 d3-contour 和 plotly.js 源码改进的等值线生成库。该实现结合了两个库的优点，并添加了额外的增强功能。

## 核心改进点

### 1. 算法基础 - Marching Squares 增强

#### 原始 d3-contour 实现
- 基础的 Marching Squares 算法
- 16 种基本情况（0-15）
- 简单的线段连接
- 基本的多边形生成

#### plotly.js 优化
- 鞍点（saddle point）消歧处理
- 更精确的等值线追踪
- 优化的内存管理
- 高效的片段拼接算法

#### enhanced-contour.js 综合改进
```javascript
// 增强的鞍点消歧函数
function disambiguateSaddle(caseIndex, corners, value) {
    if (caseIndex !== 5 && caseIndex !== 10) return caseIndex;
    
    const avg = (corners[0][0] + corners[0][1] + corners[1][0] + corners[1][1]) / 4;
    
    // 使用平均值判断鞍点类型
    if (caseIndex === 5) {
        if (value > avg) return 713; // 两个峰值带一个谷值
        return 104; // 两个谷值带一个脊线
    } else if (caseIndex === 10) {
        if (value > avg) return 1114;
        return 208;
    }
}
```

### 2. 平滑算法优化

#### d3-contour 方法
- 基本的线性插值
- 固定的平滑参数

#### enhanced-contour.js 改进
- 可配置的平滑等级（0-5）
- 智能边界检测
- 增强的插值算法
```javascript
function createSmoothingFunction(smoothingLevel = 1.0) {
    return function smoothAdvanced(ring, values, value, dx, dy) {
        const smoothingFactor = Math.max(0, Math.min(5, smoothingLevel));
        // 仅在接近边界时进行平滑
        if (Math.abs(xFrac - 0.5) < 0.1) {
            point[0] = smoothInterpolate(x, v0, v1, value, smoothingFactor);
        }
    };
}
```

### 3. 空间索引与孔洞分配

#### plotly.js 方法
- 简单的包含性测试
- O(n²) 复杂度的孔洞分配

#### enhanced-contour.js 优化
- 专用的空间索引类
- 边界框预筛选
- 两阶段孔洞分配算法
```javascript
class SpatialIndex {
    findContainingPolygon(hole) {
        // 第一阶段：边界框测试
        const candidates = this.bounds.filter(bounds => 
            holeBounds.minX >= bounds.minX && holeBounds.maxX <= bounds.maxX
        );
        
        // 第二阶段：精确包含测试
        for (const idx of candidates) {
            if (contains(this.polygons[idx].polygon[0], hole) !== -1) {
                return idx;
            }
        }
    }
}
```

### 4. 空值掩码支持

#### 新增功能
- 支持不规则数据边界
- 空值区域的特殊处理
- 自动生成边界轮廓
```javascript
contours.nullMask = function(_) {
    if (!arguments.length) return nullMask;
    nullMask = _;
    return contours;
};
```

### 5. 多模式输出

#### 支持三种输出模式
1. **lines**: 仅生成等值线（开放路径）
2. **surfaces**: 仅生成等值面（多边形）
3. **both**: 同时生成等值线和等值面

```javascript
contours.mode = function(_) {
    if (!['lines', 'surfaces', 'both'].includes(_)) 
        throw new Error("invalid mode");
    mode = _;
    return contours;
};
```

## 性能优化

### 1. 内存管理
- 使用 Map 替代普通对象进行片段索引
- 减少不必要的数组复制
- 智能的片段合并策略

### 2. 计算优化
- 预计算边界框
- 避免重复的包含性测试
- 优化的网格遍历顺序

### 3. 数据结构
- 高效的片段存储结构
- 优化的环形数据表示
- 减少临时对象创建

## GeoJSON 导出增强

### 标准 GeoJSON 支持
```javascript
export function toGeoJSON(contourResult, properties = {}) {
    // 支持单个结果或结果数组
    if (Array.isArray(contourResult)) {
        return {
            type: "FeatureCollection",
            features: contourResult.map((item, i) => 
                toGeoJSON(item, { ...properties, index: i })
            )
        };
    }
    
    // 支持混合模式输出
    if (contourResult.lines && contourResult.surfaces) {
        return {
            type: "FeatureCollection",
            features: [
                {
                    type: "Feature",
                    properties: { type: "lines", value: contourResult.value },
                    geometry: contourResult.lines
                },
                {
                    type: "Feature", 
                    properties: { type: "surfaces", value: contourResult.value },
                    geometry: contourResult.surfaces
                }
            ]
        };
    }
}
```

## 使用示例

```javascript
// 创建等值线生成器
const contour = enhancedContours()
    .size([width, height])
    .thresholds(10)  // 10个等值级别
    .smooth(2.0)     // 平滑等级2
    .mode('both')    // 同时生成线和面
    .nullMask(mask); // 设置空值掩码

// 生成等值线
const contours = contour(values);

// 导出为 GeoJSON
const geojson = toGeoJSON(contours, {
    dataset: 'testData1',
    timestamp: new Date().toISOString()
});
```

## 与原始实现的主要区别

### d3-contour
- **原始**: 基础 Marching Squares，简单平滑
- **增强**: 鞍点消歧，可配置平滑，空值支持

### plotly.js
- **原始**: 优化的片段拼接，固定输出格式
- **增强**: 空间索引优化，多模式输出，GeoJSON 导出

## 已知问题与解决方案

### 坐标为空的问题
当前实现中，如果输入数据格式不正确或阈值设置不当，可能导致生成的坐标为空。需要确保：

1. 输入数据为一维数组（按行优先顺序）
2. 正确设置 size([width, height])
3. 阈值在数据范围内

### 性能考虑
- 对于大型数据集（>10000 点），建议使用较少的阈值级别
- 可以通过降低平滑等级来提高性能
- 使用空值掩码可以跳过无效区域的计算

## 总结

`enhanced-contour.js` 成功结合了 d3-contour 的简洁性和 plotly.js 的高级特性，同时添加了额外的功能如空值掩码支持、多模式输出和 GeoJSON 导出。这使得它成为一个功能强大且灵活的等值线生成工具。