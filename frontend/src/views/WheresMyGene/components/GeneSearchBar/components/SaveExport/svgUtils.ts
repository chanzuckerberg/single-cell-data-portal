import { Tissue } from "src/views/WheresMyGene/common/types";
import {
  ECHART_AXIS_LABEL_COLOR_HEX,
  ECHART_AXIS_LABEL_FONT_SIZE_PX,
} from "../../../HeatMap/components/XAxisChart/style";
import {
  HEAT_MAP_BASE_CELL_PX,
  X_AXIS_CHART_HEIGHT_PX,
  Y_AXIS_CHART_WIDTH_PX,
} from "../../../HeatMap/utils";
import { PLASMA_SVG_STRING } from "../../../Filters/components/ColorScale";
import { CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX } from "src/components/Layout/style";
import { CHART_PADDING_PX } from "../../../HeatMap/style";
import { Y_AXIS_TISSUE_WIDTH_PX } from "../../../HeatMap/components/YAxisChart/style";
import { capitalize } from "../../../HeatMap/components/YAxisChart";

export const NAME_SPACE_URI = "http://www.w3.org/2000/svg";
const FONT_FAMILY = "sans-serif";
const LABEL_ROTATION = "rotate(-90)";

export function renderLegend({
  heatmapContainer,
  yOffset,
}: {
  heatmapContainer?: HTMLElement | null;
  yOffset: number;
}) {
  if (!heatmapContainer) return;

  const width = 120;

  const relativeGeneExpressionElement = heatmapContainer.querySelector(
    `#relative-gene-expression`
  );
  const expressedInCellsElement =
    heatmapContainer.querySelector(`#expressed-in-cells`);

  const colorScale = renderColorScale({
    element: relativeGeneExpressionElement,
    width,
  });

  const expressedInCells = renderExpressedInCells({
    element: expressedInCellsElement,
    width,
  });

  // Create root SVG element
  const svg = document.createElementNS(NAME_SPACE_URI, "svg");

  const svgAttributes = {
    fill: ECHART_AXIS_LABEL_COLOR_HEX,
    "font-family": FONT_FAMILY,
    "font-size": ECHART_AXIS_LABEL_FONT_SIZE_PX,
    x: `0`,
    y: `${yOffset}`,
  };

  applyAttributes(svg, svgAttributes);

  svg.append(colorScale || "");
  svg.append(expressedInCells || "");

  return svg;
}

export function renderExpressedInCells({
  element,
  width = 0,
}: {
  element?: Element | null;
  width?: number;
}) {
  if (!element) return;

  const expressedCellsGroup = document.createElementNS(NAME_SPACE_URI, "g");
  expressedCellsGroup.id = "expressed-in-cells";

  const xOffset = width + 60;

  // Expressed in cells label
  const expressedCellsLabel = document.createElementNS(NAME_SPACE_URI, "text");
  expressedCellsLabel.textContent =
    element.querySelector(`#expressed-in-cells-label`)?.innerHTML || null;
  applyAttributes(expressedCellsLabel, {
    transform: `translate(${xOffset}, 45)`,
  });

  // Expressed in cells dots
  const expressedCellsDots = document.createElementNS(NAME_SPACE_URI, "g");
  applyAttributes(expressedCellsDots, {
    fill: "#CCCCCC",
  });

  const dots = element?.querySelectorAll(`#expressed-in-cells-dots span`);
  Array.from(dots || []).forEach((dot, index) => {
    const circle = document.createElementNS(NAME_SPACE_URI, "circle");

    const circleAttributes = {
      r: `${Number(dot.getAttribute("size")) / 2}`,
      transform: `translate(${xOffset + 2 + 28 * index}, 60)`,
    };

    applyAttributes(circle, circleAttributes);

    expressedCellsDots.append(circle);
  });

  // Expressed in cells values
  const expressedCellsValuesGroup = document.createElementNS(
    NAME_SPACE_URI,
    "g"
  );
  applyAttributes(expressedCellsValuesGroup, {
    transform: `translate(${xOffset}, 80)`,
  });

  const expressedCellsValues = element.querySelectorAll(`.low-high span`);

  const expressedCellsValueLow = document.createElementNS(
    NAME_SPACE_URI,
    "text"
  );
  applyAttributes(expressedCellsValueLow, {
    x: 0,
  });
  expressedCellsValueLow.textContent =
    expressedCellsValues[0].innerHTML || null;

  const expressedCellsValueHigh = document.createElementNS(
    NAME_SPACE_URI,
    "text"
  );
  applyAttributes(expressedCellsValueHigh, {
    x: width,
    "text-anchor": "end",
  });
  expressedCellsValueHigh.textContent =
    expressedCellsValues[1].innerHTML || null;

  expressedCellsValuesGroup.append(expressedCellsValueLow);
  expressedCellsValuesGroup.append(expressedCellsValueHigh);

  // Append color scale legend
  expressedCellsGroup.append(expressedCellsLabel);
  expressedCellsGroup.append(expressedCellsDots);
  expressedCellsGroup.append(expressedCellsValuesGroup);

  return expressedCellsGroup;
}

export function renderColorScale({
  element,
  width = 0,
}: {
  element?: Element | null;
  width?: number;
}) {
  if (!element) return;

  const xOffset = "40";

  const colorScaleGroup = document.createElementNS(NAME_SPACE_URI, "g");
  colorScaleGroup.id = "color-scale";

  // Color scale label
  const colorScaleLabel = document.createElementNS(NAME_SPACE_URI, "text");
  colorScaleLabel.textContent =
    element.querySelector(`#relative-gene-expression-label`)?.innerHTML || null;
  applyAttributes(colorScaleLabel, {
    transform: `translate(${xOffset}, 45)`,
  });

  // Color scale image
  const colorScaleImage = document.createElementNS(NAME_SPACE_URI, "g");
  colorScaleImage.innerHTML = PLASMA_SVG_STRING;
  applyAttributes(colorScaleImage, {
    transform: `translate(${xOffset}, 50)`,
    width: width,
  });

  // Color scale values
  const colorScaleValuesGroup = document.createElementNS(NAME_SPACE_URI, "g");
  applyAttributes(colorScaleValuesGroup, {
    transform: `translate(${xOffset}, 80)`,
  });
  const colorScaleValues = element.querySelectorAll(`.low-high span`);

  const colorScaleValueLow = document.createElementNS(NAME_SPACE_URI, "text");
  applyAttributes(colorScaleValueLow, {
    x: 0,
  });
  colorScaleValueLow.textContent = colorScaleValues[0].innerHTML || null;

  const colorScaleValueHigh = document.createElementNS(NAME_SPACE_URI, "text");
  applyAttributes(colorScaleValueHigh, {
    x: 120,
    "text-anchor": "end",
  });
  colorScaleValueHigh.textContent = colorScaleValues[1].innerHTML || null;

  colorScaleValuesGroup.append(colorScaleValueLow);
  colorScaleValuesGroup.append(colorScaleValueHigh);

  // Append color scale legend
  colorScaleGroup.append(colorScaleLabel);
  colorScaleGroup.append(colorScaleImage);
  colorScaleGroup.append(colorScaleValuesGroup);

  return colorScaleGroup;
}

export function renderXAxis({
  heatmapContainer,
  yOffset,
}: {
  heatmapContainer?: HTMLElement | null;
  yOffset: number;
}) {
  if (!heatmapContainer) return;

  const xAxis = heatmapContainer.querySelector(`#x-axis-wrapper`);

  // Create root SVG element
  const svg = document.createElementNS(NAME_SPACE_URI, "svg");

  const svgAttributes = {
    id: "x-axis-container",
    fill: ECHART_AXIS_LABEL_COLOR_HEX,
    "font-family": FONT_FAMILY,
    "font-size": ECHART_AXIS_LABEL_FONT_SIZE_PX,
  };

  applyAttributes(svg, svgAttributes);

  // Create cell type count label
  const cellTypeCount = document.createElementNS(NAME_SPACE_URI, "text");

  const xOffset = Y_AXIS_CHART_WIDTH_PX + CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX;

  const cellTypeCountAttributes = {
    transform: LABEL_ROTATION,
    x: -(X_AXIS_CHART_HEIGHT_PX + yOffset + 10),
    y: xOffset,
  };

  applyAttributes(cellTypeCount, cellTypeCountAttributes);
  cellTypeCount.textContent = "Cell Count";

  svg.append(cellTypeCount);

  // Create gene labels
  const geneLabelContainer = document.createElementNS(NAME_SPACE_URI, "g");

  Array.from(
    xAxis?.querySelectorAll(`[data-testid*='gene-label-'] span`) || []
  ).forEach((label, index) => {
    const geneLabelText = document.createElementNS(NAME_SPACE_URI, "text");

    const labelAttributes = {
      transform: LABEL_ROTATION,
      "font-size": "12px",
      x: -(X_AXIS_CHART_HEIGHT_PX + yOffset + 10),
      y:
        xOffset +
        HEAT_MAP_BASE_CELL_PX +
        CHART_PADDING_PX / 2 +
        index * HEAT_MAP_BASE_CELL_PX,
    };

    applyAttributes(geneLabelText, labelAttributes);
    geneLabelText.textContent = label.innerHTML;

    geneLabelContainer.append(geneLabelText);
  });

  svg.append(geneLabelContainer);

  return svg;
}

export function renderYAxis({
  heatmapContainer,
  tissueName,
  heatmapHeight,
  yOffset,
}: {
  heatmapContainer?: HTMLElement | null;
  tissueName: Tissue;
  heatmapHeight: number;
  yOffset: number;
}) {
  if (!heatmapContainer) return;

  const yAxis = heatmapContainer.querySelector(`#${tissueName}-y-axis`);

  const svg = document.createElementNS(NAME_SPACE_URI, "svg");

  const svgAttributes = {
    height: `${heatmapHeight}px`,
    width: `${Y_AXIS_CHART_WIDTH_PX}px`, // Adds padding for current tissue label
    x: `${CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX}`,
    y: `${X_AXIS_CHART_HEIGHT_PX + yOffset + 25}`,
  };

  applyAttributes(svg, svgAttributes);

  // Container for tissue label
  const tissueLabelContainer = document.createElementNS(NAME_SPACE_URI, "svg");

  const RECT_WIDTH = 5;

  const tissueLabelText = document.createElementNS(NAME_SPACE_URI, "text");
  tissueLabelText.textContent = capitalize(tissueName);
  applyAttributes(tissueLabelText, {
    y: `${Y_AXIS_TISSUE_WIDTH_PX - RECT_WIDTH - 12}`,
    "text-anchor": "end",
    transform: LABEL_ROTATION,
    "font-family": "Inter, sans-serif",
    "font-size": `14px`,
    "font-weight": "bold",
  });

  const tissueLabelRect = document.createElementNS(NAME_SPACE_URI, "rect");
  applyAttributes(tissueLabelRect, {
    width: `${RECT_WIDTH}`,
    height: `${heatmapHeight}`,
    x: `${Y_AXIS_TISSUE_WIDTH_PX - RECT_WIDTH}`,
    fill: "black",
  });

  tissueLabelContainer.append(tissueLabelText);
  tissueLabelContainer.append(tissueLabelRect);

  // cell type names container group
  const cellTypeNamesContainer = document.createElementNS(
    NAME_SPACE_URI,
    "svg"
  );
  applyAttributes(cellTypeNamesContainer, {
    x: `${Y_AXIS_TISSUE_WIDTH_PX + 10}`,
  });

  Array.from(
    yAxis?.querySelectorAll("[data-testid='cell-type-label-count']") || []
  ).forEach((labelCount, index) => {
    const label = labelCount.querySelector(
      "[data-testid='cell-type-name']"
    )?.textContent;

    const count = labelCount.querySelector(
      "[data-testid='cell-count']"
    )?.textContent;

    // group
    const group = document.createElementNS(NAME_SPACE_URI, "svg");

    const groupAttributes = {
      id: "cell-name-label-group",
      fill: ECHART_AXIS_LABEL_COLOR_HEX,
      "font-family": FONT_FAMILY,
      "font-size": ECHART_AXIS_LABEL_FONT_SIZE_PX,
      x: 0,
      /**
       * (thuang): Add `HEAT_MAP_BASE_CELL_PX / 2` top margin, so we render the
       * first label in the middle of the first cell
       */
      y: index * HEAT_MAP_BASE_CELL_PX,
    };

    applyAttributes(group, groupAttributes);

    // cellTypeLabel
    const cellTypeLabel = document.createElementNS(NAME_SPACE_URI, "text");

    // Preserves whitespace
    cellTypeLabel.setAttributeNS(
      "http://www.w3.org/XML/1998/namespace",
      "xml:space",
      "preserve"
    );

    const labelAttributes = {
      x: 0,
      y: HEAT_MAP_BASE_CELL_PX / 2,
    };

    applyAttributes(cellTypeLabel, labelAttributes);
    cellTypeLabel.textContent = String(label);

    // cellCount
    const cellCount = document.createElementNS(NAME_SPACE_URI, "text");

    const countAttributes = {
      "text-anchor": "end",
      x: Y_AXIS_CHART_WIDTH_PX - CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX,
      y: HEAT_MAP_BASE_CELL_PX / 2,
    };

    applyAttributes(cellCount, countAttributes);

    cellCount.textContent = String(count);

    // append children
    group.appendChild(cellTypeLabel);
    group.appendChild(cellCount);

    cellTypeNamesContainer.appendChild(group);
  });

  svg.append(tissueLabelContainer);
  svg.append(cellTypeNamesContainer);

  return svg;
}

export function renderDots({
  heatmapContainer,
  tissueName,
  yOffset,
}: {
  heatmapContainer?: HTMLElement | null;
  tissueName: Tissue;
  yOffset: number;
}) {
  if (!heatmapContainer) return;

  const chart = heatmapContainer
    .querySelector(`#${tissueName}-chart`)
    ?.querySelector("svg");

  if (!chart) return;

  applyAttributes(chart, {
    x:
      Y_AXIS_CHART_WIDTH_PX +
      CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX +
      CHART_PADDING_PX,
    y: X_AXIS_CHART_HEIGHT_PX + yOffset + 20,
  });

  // Cleanup as style attributes aren't used in SVG files
  chart.removeAttribute("style");
  Array.from(chart?.querySelectorAll("rect, path, g") || []).forEach(
    (element) => {
      element.removeAttribute("style");
    }
  );

  return chart;
}

export function applyAttributes(
  node: SVGElement | HTMLElement,
  attributes: Record<string, string | number>
) {
  for (const [key, value] of Object.entries(attributes)) {
    node.setAttribute(key, String(value));
  }
}
