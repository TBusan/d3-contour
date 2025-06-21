/**
 * 颜色比例尺工具
 * 用于生成颜色映射
 */

// 预定义的颜色方案
const COLOR_SCHEMES = {
  // 蓝色到红色
  'blueRed': [
    [0, 0, 255],    // 蓝色
    [255, 0, 0]     // 红色
  ],
  
  // 彩虹色
  'rainbow': [
    [148, 0, 211],  // 紫色
    [75, 0, 130],   // 靛色
    [0, 0, 255],    // 蓝色
    [0, 255, 0],    // 绿色
    [255, 255, 0],  // 黄色
    [255, 127, 0],  // 橙色
    [255, 0, 0]     // 红色
  ],
  
  // Viridis (来自Matplotlib)
  'viridis': [
    [68, 1, 84],
    [65, 68, 135],
    [42, 120, 142],
    [34, 168, 132],
    [122, 209, 81],
    [253, 231, 37]
  ],
  
  // Plasma (来自Matplotlib)
  'plasma': [
    [13, 8, 135],
    [126, 3, 168],
    [204, 71, 120],
    [248, 149, 64],
    [240, 249, 33]
  ],
  
  // Inferno (来自Matplotlib)
  'inferno': [
    [0, 0, 4],
    [51, 13, 53],
    [122, 28, 69],
    [192, 82, 64],
    [252, 180, 45],
    [252, 253, 191]
  ],
  
  // Magma (来自Matplotlib)
  'magma': [
    [0, 0, 4],
    [43, 13, 53],
    [120, 28, 109],
    [187, 55, 84],
    [249, 142, 82],
    [252, 253, 191]
  ],
  
  // 蓝色到白色
  'blueWhite': [
    [0, 0, 255],    // 蓝色
    [255, 255, 255] // 白色
  ],
  
  // 热力图
  'hot': [
    [0, 0, 0],      // 黑色
    [255, 0, 0],    // 红色
    [255, 255, 0],  // 黄色
    [255, 255, 255] // 白色
  ],
  
  // 冷色调
  'cool': [
    [0, 255, 255],  // 青色
    [255, 0, 255]   // 洋红色
  ],
  
  // 蓝绿红
  'bgr': [
    [0, 0, 255],    // 蓝色
    [0, 255, 0],    // 绿色
    [255, 0, 0]     // 红色
  ]
};

/**
 * 创建颜色查找表
 * @param {String|Array|Function} colorScale 颜色比例尺
 * @returns {Function} 颜色查找函数
 */
export function createLookupTable(colorScale) {
  // 如果已经是函数，直接返回
  if (typeof colorScale === 'function') {
    return colorScale;
  }
  
  // 如果是字符串，查找预定义的颜色方案
  if (typeof colorScale === 'string') {
    const scheme = COLOR_SCHEMES[colorScale];
    if (scheme) {
      return createColorInterpolator(scheme);
    }
    // 如果找不到预定义方案，使用默认方案
    return createColorInterpolator(COLOR_SCHEMES.viridis);
  }
  
  // 如果是数组，创建自定义颜色插值器
  if (Array.isArray(colorScale)) {
    return createColorInterpolator(colorScale);
  }
  
  // 默认返回viridis方案
  return createColorInterpolator(COLOR_SCHEMES.viridis);
}

/**
 * 创建颜色插值器
 * @param {Array<Array<Number>>} colors 颜色数组
 * @returns {Function} 颜色插值函数
 */
function createColorInterpolator(colors) {
  return function(t) {
    // 确保t在[0,1]范围内
    t = Math.max(0, Math.min(1, t));
    
    // 如果只有一种颜色，直接返回
    if (colors.length === 1) {
      return rgbToString(colors[0]);
    }
    
    // 计算t在颜色数组中的位置
    const index = t * (colors.length - 1);
    const i = Math.floor(index);
    const f = index - i; // 小数部分
    
    // 如果t刚好落在某个颜色上，直接返回该颜色
    if (f === 0) {
      return rgbToString(colors[i]);
    }
    
    // 在两个颜色之间插值
    const color1 = colors[i];
    const color2 = colors[i + 1];
    
    const r = Math.round(color1[0] * (1 - f) + color2[0] * f);
    const g = Math.round(color1[1] * (1 - f) + color2[1] * f);
    const b = Math.round(color1[2] * (1 - f) + color2[2] * f);
    
    return rgbToString([r, g, b]);
  };
}

/**
 * 将RGB数组转换为CSS颜色字符串
 * @param {Array<Number>} rgb RGB数组
 * @returns {String} CSS颜色字符串
 */
function rgbToString(rgb) {
  return `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
}

/**
 * 创建带透明度的颜色
 * @param {String} color 颜色字符串
 * @param {Number} opacity 透明度（0-1）
 * @returns {String} 带透明度的颜色字符串
 */
export function createColorWithOpacity(color, opacity) {
  // 如果颜色已经是rgba格式，替换透明度
  if (color.startsWith('rgba')) {
    return color.replace(/rgba\((.+?), .+?\)/, `rgba($1, ${opacity})`);
  }
  
  // 如果颜色是rgb格式，转换为rgba
  if (color.startsWith('rgb')) {
    return color.replace(/rgb\((.+?)\)/, `rgba($1, ${opacity})`);
  }
  
  // 如果颜色是十六进制格式，转换为rgba
  if (color.startsWith('#')) {
    const r = parseInt(color.slice(1, 3), 16);
    const g = parseInt(color.slice(3, 5), 16);
    const b = parseInt(color.slice(5, 7), 16);
    return `rgba(${r}, ${g}, ${b}, ${opacity})`;
  }
  
  // 默认返回带透明度的黑色
  return `rgba(0, 0, 0, ${opacity})`;
}

/**
 * 获取可用的颜色方案名称
 * @returns {Array<String>} 颜色方案名称数组
 */
export function getAvailableColorSchemes() {
  return Object.keys(COLOR_SCHEMES);
} 