import assert from "assert";
import {
  padGrid,
  chaikinSmooth,
  catmullRomSmooth,
  autoLevels,
  generateContours,
  generateContoursAuto,
  extractIsolinesFromBands,
  generateLabelsForContours,
  removeCollidingLabels
} from "../src/index.js";

describe("Enhanced d3-contour", () => {
  describe("padGrid", () => {
    it("should pad a 2x2 grid to 4x4", () => {
      const grid = new Float32Array([1, 2, 3, 4]);
      const padded = padGrid(grid, 2, 2);

      assert.equal(padded.length, 16);
      assert.equal(padded[1 * 4 + 1], 1); // Original value
      assert.equal(padded[1 * 4 + 2], 2);
      assert.equal(padded[2 * 4 + 1], 3);
      assert.equal(padded[2 * 4 + 2], 4);
    });

    it("should pad edges correctly", () => {
      const grid = new Float32Array([1, 2, 3, 4]);
      const padded = padGrid(grid, 2, 2);

      // Top edge should match row 1
      assert.equal(padded[0 * 4 + 1], 1);
      assert.equal(padded[0 * 4 + 2], 2);
    });
  });

  describe("chaikinSmooth", () => {
    it("should smooth a simple path", () => {
      const points = [[0, 0], [1, 1], [2, 0]];
      const smoothed = chaikinSmooth(points, 1);

      assert.ok(smoothed.length > points.length);
    });

    it("should handle closed paths", () => {
      const points = [[0, 0], [1, 1], [2, 0], [0, 0]];
      const smoothed = chaikinSmooth(points, 1);

      assert.ok(smoothed.length > points.length);
    });
  });

  describe("catmullRomSmooth", () => {
    it("should smooth a path", () => {
      const points = [[0, 0], [1, 0.5], [2, 0], [3, 0.5], [4, 0]];
      const smoothed = catmullRomSmooth(points, 0.5, 10);

      assert.ok(smoothed.length >= points.length);
    });
  });

  describe("autoLevels", () => {
    it("should generate nice levels", () => {
      const levels = autoLevels(3.7, 97.3, { count: 10 });

      assert.ok(levels.length > 0);
      assert.ok(levels[0] >= 0);
      assert.ok(levels[levels.length - 1] <= 100);
    });

    it("should generate round numbers", () => {
      const levels = autoLevels(0, 100, { count: 5 });

      // Check that levels are multiples of 10, 20, 50, or 100
      const step = levels[1] - levels[0];
      assert.ok([10, 20, 25, 50].includes(step));
    });
  });

  describe("generateContours", () => {
    it("should generate bands and lines", () => {
      const nx = 10, ny = 10;
      const grid = new Float32Array(nx * ny);

      for (let y = 0; y < ny; y++) {
        for (let x = 0; x < nx; x++) {
          grid[y * nx + x] = Math.sin(x / 3) * Math.cos(y / 3);
        }
      }

      const { bands, lines } = generateContours(grid, {
        size: [nx, ny],
        thresholds: [-0.5, 0, 0.5]
      });

      assert.ok(Array.isArray(bands));
      assert.ok(Array.isArray(lines));
    });

    it("should generate bands with correct structure", () => {
      const nx = 10, ny = 10;
      const grid = new Float32Array(nx * ny);

      for (let i = 0; i < grid.length; i++) {
        grid[i] = (i % 5) / 5;
      }

      const { bands } = generateContours(grid, {
        size: [nx, ny],
        thresholds: [0.2, 0.4, 0.6]
      });

      assert.ok(bands.length > 0);
      assert.ok(bands[0].value); // [min, max] range
      assert.ok(bands[0].polygons);
      assert.ok(typeof bands[0].level === 'number');
    });
  });

  describe("generateContoursAuto", () => {
    it("should generate contours with automatic levels", () => {
      const nx = 10, ny = 10;
      const grid = new Float32Array(nx * ny);

      for (let y = 0; y < ny; y++) {
        for (let x = 0; x < nx; x++) {
          grid[y * nx + x] = Math.random();
        }
      }

      const { bands, lines, thresholds } = generateContoursAuto(grid, {
        size: [nx, ny],
        count: 5
      });

      assert.ok(Array.isArray(bands));
      assert.ok(Array.isArray(lines));
      assert.ok(Array.isArray(thresholds));
      assert.equal(thresholds.length, 5);
    });
  });

  describe("extractIsolinesFromBands", () => {
    it("should extract isolines from bands", () => {
      const nx = 10, ny = 10;
      const grid = new Float32Array(nx * ny);

      for (let i = 0; i < grid.length; i++) {
        grid[i] = Math.random();
      }

      const { bands } = generateContours(grid, {
        size: [nx, ny],
        thresholds: [0.3, 0.6]
      });

      const lines = extractIsolinesFromBands(bands);

      assert.ok(Array.isArray(lines));
      if (lines.length > 0) {
        assert.ok(lines[0].paths);
        assert.ok(typeof lines[0].level === 'number');
      }
    });
  });

  describe("generateLabelsForContours", () => {
    it("should generate labels for contours", () => {
      const nx = 20, ny = 20;
      const grid = new Float32Array(nx * ny);

      for (let y = 0; y < ny; y++) {
        for (let x = 0; x < nx; x++) {
          grid[y * nx + x] = Math.sin(x / 5) * Math.cos(y / 5);
        }
      }

      const { lines } = generateContours(grid, {
        size: [nx, ny],
        thresholds: [0.2, 0.4, 0.6]
      });

      const labels = generateLabelsForContours(lines);

      assert.ok(Array.isArray(labels));
    });

    it("should filter colliding labels", () => {
      const labels = [
        { x: 0, y: 0, angle: 0, text: '1.0', level: 1 },
        { x: 10, y: 0, angle: 0, text: '2.0', level: 2 },
        { x: 5, y: 0, angle: 0, text: '1.5', level: 1.5 } // Collides with first two
      ];

      const filtered = removeCollidingLabels(labels, { width: 20, height: 10 });

      assert.ok(filtered.length < labels.length);
    });
  });

  describe("integration test", () => {
    it("should work end-to-end with padding, smoothing, and labels", () => {
      const nx = 30, ny = 30;
      const grid = new Float32Array(nx * ny);

      // Create a simple gaussian
      for (let y = 0; y < ny; y++) {
        for (let x = 0; x < nx; x++) {
          const dx = x - nx / 2;
          const dy = y - ny / 2;
          grid[y * nx + x] = Math.exp(-(dx * dx + dy * dy) / 50);
        }
      }

      // Pad grid
      const padded = padGrid(grid, nx, ny);

      // Generate contours
      const { bands, lines } = generateContoursAuto(padded, {
        size: [nx + 2, ny + 2],
        count: 5
      });

      // Smooth lines
      const smoothedLines = lines.map(line => ({
        ...line,
        paths: line.paths.map(path => chaikinSmooth(path, 2))
      }));

      // Generate labels
      const labels = removeCollidingLabels(
        generateLabelsForContours(smoothedLines),
        { width: 30, height: 15 }
      );

      // Verify results
      assert.ok(bands.length > 0);
      assert.ok(lines.length > 0);
      assert.ok(Array.isArray(labels));
    });
  });
});
