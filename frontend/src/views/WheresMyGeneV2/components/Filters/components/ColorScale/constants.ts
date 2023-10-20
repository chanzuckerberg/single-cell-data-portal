// (ashin-czi): Used by SaveExport SVG generation to recreate the color scale plasma image in an SVG
// Not the best solution but importing SVG files creates conflicts with the Loader

import { SORT_BY } from "src/views/WheresMyGene/common/types";

// and can't get the actual content of the SVG file needed for SVG creation
export const PLASMA_SVG_STRING = `
  <rect width="120" height="16" fill="url(#paint0_linear_5151_537460)" />
  <defs>
      <linearGradient id="paint0_linear_5151_537460" x1="120" y1="8" x2="8.0516e-07" y2="8.00001"
          gradientUnits="userSpaceOnUse">
          <stop stop-color="#090720" />
          <stop offset="0.145123" stop-color="#36106B" />
          <stop offset="0.28796" stop-color="#6B1D81" />
          <stop offset="0.436206" stop-color="#9C2E7F" />
          <stop offset="0.582062" stop-color="#D3436E" />
          <stop offset="0.719889" stop-color="#F66E5C" />
          <stop offset="0.864407" stop-color="#FEA973" />
          <stop offset="1" stop-color="#FDE2A3" />
      </linearGradient>
  </defs>
`;

export const COLOR_SCALE_OPTIONS = [
  { id: SORT_BY.COLOR_SCALED, name: "Scaled" },
  { id: SORT_BY.COLOR_UNSCALED, name: "Unscaled" },
];
