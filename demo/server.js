/**
 * Simple HTTP server for serving the demo
 * Run with: node demo/server.js
 */

import { createServer } from 'http';
import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

const PORT = 8080;

const mimeTypes = {
    '.html': 'text/html',
    '.js': 'application/javascript',
    '.css': 'text/css',
    '.json': 'application/json',
    '.png': 'image/png',
    '.jpg': 'image/jpg',
    '.gif': 'image/gif',
    '.svg': 'image/svg+xml',
    '.ico': 'image/x-icon'
};

const server = createServer((req, res) => {
    console.log(`${req.method} ${req.url}`);

    // Default to index.html
    let filePath = '.' + req.url;
    if (filePath === './') {
        filePath = './demo/interactive-contour.html';
    }

    // Get file extension
    const extname = String(filePath).split('.').pop().toLowerCase();
    const contentType = mimeTypes['.' + extname] || 'application/octet-stream';

    // Read and serve file
    try {
        const content = readFileSync(filePath);
        res.writeHead(200, {
            'Content-Type': contentType,
            'Access-Control-Allow-Origin': '*'
        });
        res.end(content, 'utf-8');
    } catch (err) {
        if (err.code === 'ENOENT') {
            res.writeHead(404, { 'Content-Type': 'text/html' });
            res.end('<h1>404 Not Found</h1>', 'utf-8');
        } else {
            res.writeHead(500);
            res.end(`Server Error: ${err.code}`, 'utf-8');
        }
    }
});

server.listen(PORT, () => {
    console.log('\n=================================');
    console.log(`🚀 Demo server running at:`);
    console.log(`   http://localhost:${PORT}`);
    console.log('=================================\n');
    console.log('Press Ctrl+C to stop the server\n');
});
