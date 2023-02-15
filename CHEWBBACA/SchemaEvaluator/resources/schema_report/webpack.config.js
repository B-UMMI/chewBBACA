var path = require("path");

// webpack needs to be explicitly required
const webpack = require("webpack");

module.exports = {
  mode: "development",
  entry: "./src/index.js",
  output: {
    path: path.resolve(__dirname, "dist"),
    filename: "main.js",
  },
  resolve: {
    fallback: { 'process/browser': require.resolve('process/browser'), }
  },
  module: {
    rules: [
      { test: /\.(js|jsx)$/, exclude: /node_modules/, use: "babel-loader" },
      {
        test: /\.scss$/,
        use: ["style-loader", "css-loader", "sass-loader"],
      },
      {
        test: /\.css$/,
        use: [
          "style-loader",
          {
            loader: "css-loader",
            options: {
              importLoaders: 1,
              modules: true,
            },
          },
        ],
      },
      {
        test: /\.(png|jpe?g|gif)$/i,
        use: [
          {
            loader: "file-loader",
          },
        ],
      },
      {
        test: /\.md$/,
        use: [
          {
            loader: "remark-loader",
          },
        ],
      },
      {
        test: /\.(fasta|clustal|nh)$/i,
        use: 'raw-loader',
      },
    ],
  },
  plugins: [
    // fix "process is not defined" error:
    new webpack.ProvidePlugin({
      process: "process/browser",
    }),
  ],
};
