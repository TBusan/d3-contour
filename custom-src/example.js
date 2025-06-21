import { contours, contourBands } from "./index.js";

// Example function to demonstrate the usage of our custom contour library
function demo() {
  // Example data with null values (5x5 grid)
  const data = [
    [0.5, 1.0, 1.5, 2.0, 1.5],
    [1.0, null, 2.5, 3.0, 2.0],
    [1.5, 2.5, 4.0, null, 2.5],
    [2.0, 3.0, 3.5, 3.0, 2.0],
    [1.5, 2.0, 2.5, 2.0, 1.5]
  ];
  
  // Flatten the 2D array for the contour generator
  const flatData = [];
  for (const row of data) {
    for (const value of row) {
      flatData.push(value);
    }
  }
  
  console.log("Generating contours with null values...");

  // Create and configure contour generator
  const contourGenerator = contours()
    .size([5, 5]) // Width x height
    .smooth(true)
    .thresholds([1, 2, 3, 4]); // Custom threshold values
  
  // Generate contours
  const contourData = contourGenerator(flatData);
  console.log(`Generated ${contourData.length} contours`);
  
  // Log the contour values
  contourData.forEach(c => {
    console.log(`Contour at value ${c.value}: ${c.coordinates.length} polygon(s)`);
  });

  console.log("\nGenerating contour bands...");
  
  // Create and configure contour band generator
  const bandGenerator = contourBands()
    .size([5, 5])
    .smooth(true)
    .thresholds([1, 2, 3, 4, 5]);
  
  // Generate bands
  const bandData = bandGenerator(flatData);
  console.log(`Generated ${bandData.length} bands`);
  
  // Log the band values
  bandData.forEach(b => {
    console.log(`Band from ${b.lowerValue} to ${b.upperValue}: ${b.coordinates.length} polygon(s)`);
  });
  
  return {
    contours: contourData,
    bands: bandData
  };
}

// Run the demo if this file is executed directly
if (typeof window !== 'undefined') {
  window.runDemo = demo;
  console.log("Demo function available as window.runDemo()");
} else if (typeof require !== 'undefined') {
  demo();
}

export { demo }; 