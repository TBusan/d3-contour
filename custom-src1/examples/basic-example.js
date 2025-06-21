/**
 * 基本示例
 * 展示如何使用等值线和等值面渲染
 */

import { generateContours, generateContourBands, generateFilledContours } from '../index.js';

// 创建一个示例数据集
function createSampleData(width, height) {
  const data = new Array(height);
  
  for (let y = 0; y < height; y++) {
    data[y] = new Array(width);
    for (let x = 0; x < width; x++) {
      // 创建一个简单的高斯分布
      const dx = x - width / 2;
      const dy = y - height / 2;
      const distance = Math.sqrt(dx * dx + dy * dy);
      
      // 添加一些波纹效果
      const wave = Math.sin(distance / 5) * 10;
      
      // 最终值
      data[y][x] = 100 - distance + wave;
    }
  }
  
  return data;
}

// 渲染等值线
function renderContourLines(data, thresholds, options = {}) {
  console.log('Generating contour lines...');
  
  // 生成等值线
  const contours = generateContours(data, thresholds, options);
  
  console.log(`Generated ${Array.isArray(contours) ? contours.length : 1} contour line(s)`);
  
  // 在这里可以使用SVG或Canvas渲染等值线
  // 例如，使用D3.js渲染到SVG：
  /*
  const svg = d3.select('#contour-lines')
    .append('svg')
    .attr('width', width)
    .attr('height', height);
    
  const contourGroup = svg.append('g');
  
  contourGroup.selectAll('path')
    .data(Array.isArray(contours) ? contours : [contours])
    .enter()
    .append('path')
    .attr('d', d => d3.geoPath()(d))
    .attr('fill', 'none')
    .attr('stroke', 'black')
    .attr('stroke-width', 1);
  */
  
  return contours;
}

// 渲染等值面
function renderContourBands(data, thresholds, options = {}) {
  console.log('Generating contour bands...');
  
  // 生成等值面
  const bands = generateContourBands(data, thresholds, options);
  
  console.log(`Generated ${Array.isArray(bands) ? bands.length : 1} contour band(s)`);
  
  // 在这里可以使用SVG或Canvas渲染等值面
  // 例如，使用D3.js渲染到SVG：
  /*
  const svg = d3.select('#contour-bands')
    .append('svg')
    .attr('width', width)
    .attr('height', height);
    
  const bandsGroup = svg.append('g');
  
  bandsGroup.selectAll('path')
    .data(Array.isArray(bands) ? bands : [bands])
    .enter()
    .append('path')
    .attr('d', d => d3.geoPath()(d))
    .attr('fill', d => d.fill.color)
    .attr('fill-opacity', d => d.fill.opacity)
    .attr('stroke', d => d.stroke ? d.stroke.color : 'none')
    .attr('stroke-width', d => d.stroke ? d.stroke.width : 0);
  */
  
  return bands;
}

// 渲染填充等值线
function renderFilledContours(data, thresholds, options = {}) {
  console.log('Generating filled contours...');
  
  // 生成填充等值线
  const filled = generateFilledContours(data, thresholds, options);
  
  console.log('Generated filled contours');
  
  // 在这里可以使用SVG或Canvas渲染填充等值线
  // 例如，使用D3.js渲染到SVG：
  /*
  const svg = d3.select('#filled-contours')
    .append('svg')
    .attr('width', width)
    .attr('height', height);
    
  const filledGroup = svg.append('g');
  
  // 渲染填充区域
  filledGroup.selectAll('path.fill')
    .data(filled.fills)
    .enter()
    .append('path')
    .attr('class', 'fill')
    .attr('d', (d, i) => d3.geoPath()(
      { type: 'MultiPolygon', coordinates: [filled.coordinates[i]] }
    ))
    .attr('fill', d => d.color)
    .attr('fill-opacity', d => d.opacity);
    
  // 渲染轮廓线
  if (filled.lines) {
    filledGroup.selectAll('path.line')
      .data(filled.lines)
      .enter()
      .append('path')
      .attr('class', 'line')
      .attr('d', d => d3.geoPath()({ type: 'MultiLineString', coordinates: d.coordinates }))
      .attr('fill', 'none')
      .attr('stroke', d => d.color)
      .attr('stroke-width', d => d.width);
  }
  */
  
  return filled;
}

// 主函数
function main() {
  console.log('Running contour example...');
  
  // 创建示例数据
  const width = 100;
  const height = 100;
  const data = createSampleData(width, height);
  
  // 定义阈值
  const thresholds = [10, 20, 30, 40, 50, 60, 70, 80, 90];
  
  // 渲染等值线
  const contours = renderContourLines(data, thresholds, {
    smooth: true,
    smoothFactor: 0.2
  });
  
  // 渲染等值面
  const bands = renderContourBands(data, thresholds, {
    colorScale: 'viridis',
    fillOpacity: 0.6,
    showLines: true
  });
  
  // 渲染填充等值线
  const filled = renderFilledContours(data, thresholds, {
    colorScale: 'plasma',
    fillMode: 'tonext',
    showLines: true
  });
  
  console.log('Example completed');
  
  return {
    data,
    contours,
    bands,
    filled
  };
}

// 运行示例
main(); 