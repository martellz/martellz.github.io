const path = require('path');

module.exports = {
  mode: 'development',
  entry: './js/MoCap.js',
  devtool: 'inline-source-map',
  output: {
    filename: 'bundle.js',
    path: path.resolve(__dirname, 'dist')
  }
};
