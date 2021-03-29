const babelOptions = {
  presets: ["@babel/preset-react", "@babel/preset-typescript"],
};

module.exports = require("babel-jest").createTransformer(babelOptions);
