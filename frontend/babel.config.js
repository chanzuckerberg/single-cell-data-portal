module.exports = {
  plugins: [
    [
      "styled-components",
      {
        displayName: true,
        preprocess: false,
        ssr: true,
      },
    ],
  ],
  presets: ["next/babel"],
};
