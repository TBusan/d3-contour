import enhancedContours, { toGeoJSON } from './enhanced-contour.js';
import { contourData } from './testData1.js';
import fs from 'fs';

console.log('Generating contours from testData1.js...\n');

// Process TestData1
const testData = {
    width: contourData.data.x.length,
    height: contourData.data.y.length,
    x: contourData.data.x,
    y: contourData.data.y,
    values: []
};

console.log('Input data structure:', {
    xLength: contourData.data.x.length,
    yLength: contourData.data.y.length,
    vRows: contourData.data.v.length,
    vCols: contourData.data.v[0] ? contourData.data.v[0].length : 0
});

// Convert 2D array to flat 1D array (rows first)
for (let j = 0; j < testData.height; j++) {
    for (let i = 0; i < testData.width; i++) {
        const value = contourData.data.v[j][i];
        testData.values.push(value === null ? NaN : value);
    }
}

// Create null mask
const nullMask = new Array(testData.values.length);
for (let i = 0; i < testData.values.length; i++) {
    nullMask[i] = isNaN(testData.values[i]) ? 0 : 1;
}

// Check data statistics
const validValues = testData.values.filter(v => !isNaN(v));
const minValue = Math.min(...validValues);
const maxValue = Math.max(...validValues);

console.log(`Data dimensions: ${testData.width}x${testData.height}`);
console.log(`Valid values: ${validValues.length}/${testData.values.length}`);
console.log(`Value range: ${minValue.toFixed(2)} to ${maxValue.toFixed(2)} (actual)`);
console.log(`Value range: ${contourData.data.zmin} to ${contourData.data.zmax} (metadata)`);

// Calculate better thresholds based on actual data
const thresholds = [];
for (let i = 0; i < 10; i++) {
    thresholds.push(minValue + (maxValue - minValue) * (i + 1) / 11);
}

// Generate contour lines
console.log('\nGenerating contour lines...');
const lineContour = enhancedContours()
    .size([testData.width, testData.height])
    .thresholds(thresholds)
    .smooth(1.0)
    .mode('lines')
    .nullMask(nullMask);

const contourLines = lineContour(testData.values);
const linesGeoJSON = toGeoJSON(contourLines, { 
    dataset: 'testData1',
    type: 'lines',
    zmin: contourData.data.zmin,
    zmax: contourData.data.zmax,
    timestamp: new Date().toISOString()
});

fs.writeFileSync('testData1-lines.geojson', JSON.stringify(linesGeoJSON, null, 2));
console.log(`✓ Saved testData1-lines.geojson (${contourLines.length} levels)`);

// Generate contour surfaces
console.log('\nGenerating contour surfaces...');
const surfaceContour = enhancedContours()
    .size([testData.width, testData.height])
    .thresholds(thresholds)
    .smooth(1.0)
    .mode('surfaces')
    .nullMask(nullMask);

const contourSurfaces = surfaceContour(testData.values);
const surfacesGeoJSON = toGeoJSON(contourSurfaces, {
    dataset: 'testData1',
    type: 'surfaces',
    zmin: contourData.data.zmin,
    zmax: contourData.data.zmax,
    timestamp: new Date().toISOString()
});

fs.writeFileSync('testData1-surfaces.geojson', JSON.stringify(surfacesGeoJSON, null, 2));
console.log(`✓ Saved testData1-surfaces.geojson (${contourSurfaces.length} levels)`);

// Generate both lines and surfaces
console.log('\nGenerating both lines and surfaces...');
const bothContour = enhancedContours()
    .size([testData.width, testData.height])
    .thresholds(thresholds)
    .smooth(1.0)
    .mode('both')
    .nullMask(nullMask);

const contourBoth = bothContour(testData.values);
const bothGeoJSON = toGeoJSON(contourBoth, {
    dataset: 'testData1',
    type: 'both',
    zmin: contourData.data.zmin,
    zmax: contourData.data.zmax,
    timestamp: new Date().toISOString()
});

fs.writeFileSync('testData1-both.geojson', JSON.stringify(bothGeoJSON, null, 2));
console.log(`✓ Saved testData1-both.geojson (${contourBoth.length} levels)`);

// Summary report
console.log('\nSummary:');
contourLines.forEach((level, i) => {
    const lineCount = level.coordinates.length;
    console.log(`  Level ${i+1}: value=${level.value.toFixed(2)}, lines=${lineCount}`);
});

console.log('\nAll GeoJSON files generated successfully!');