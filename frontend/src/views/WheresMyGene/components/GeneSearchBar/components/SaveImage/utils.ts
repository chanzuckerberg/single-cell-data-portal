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
import { plasmaSvgString } from "../../../InfoPanel/components/ColorScale";

export const NAME_SPACE_URI = "http://www.w3.org/2000/svg";
const FONT_FAMILY = "sans-serif";

export function renderLegend({
  heatmapContainer,
}: {
  heatmapContainer?: HTMLElement | null;
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
    y: `0`,
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

  const xPosition = width + 60;

  // Expressed in cells label
  const expressedCellsLabel = document.createElementNS(NAME_SPACE_URI, "text");
  expressedCellsLabel.textContent =
    element.querySelector(`#expressed-in-cells-label`)?.innerHTML || null;
  expressedCellsLabel.setAttribute("transform", `translate(${xPosition}, 45)`);

  // Expressed in cells dots
  const expressedCellsDots = document.createElementNS(NAME_SPACE_URI, "g");
  expressedCellsDots.setAttribute("fill", "#CCCCCC");

  const dots = element?.querySelectorAll(`#expressed-in-cells-dots span`);
  Array.from(dots || []).forEach((dot, index) => {
    const circle = document.createElementNS(NAME_SPACE_URI, "circle");

    const circleAttributes = {
      r: `${Number(dot.getAttribute("size")) / 2}`,
      transform: `translate(${xPosition + 2 + 28 * index}, 60)`,
    };

    applyAttributes(circle, circleAttributes);

    expressedCellsDots.append(circle);
  });

  // Expressed in cells values
  const expressedCellsValuesGroup = document.createElementNS(
    NAME_SPACE_URI,
    "g"
  );
  expressedCellsValuesGroup.setAttribute(
    "transform",
    `translate(${xPosition}, 80)`
  );
  const expressedCellsValues = element.querySelectorAll(`.low-high span`);

  const expressedCellsValueLow = document.createElementNS(
    NAME_SPACE_URI,
    "text"
  );
  expressedCellsValueLow.setAttribute("x", "0");
  expressedCellsValueLow.textContent =
    expressedCellsValues[0].innerHTML || null;

  const expressedCellsValueHigh = document.createElementNS(
    NAME_SPACE_URI,
    "text"
  );
  expressedCellsValueHigh.setAttribute("x", `${width}`);
  expressedCellsValueHigh.setAttribute("text-anchor", "end");
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

  const xPosition = "40";

  const colorScaleGroup = document.createElementNS(NAME_SPACE_URI, "g");
  colorScaleGroup.id = "color-scale";

  // Color scale label
  const colorScaleLabel = document.createElementNS(NAME_SPACE_URI, "text");
  colorScaleLabel.textContent =
    element.querySelector(`#relative-gene-expression-label`)?.innerHTML || null;
  colorScaleLabel.setAttribute("transform", `translate(${xPosition}, 45)`);

  // Color scale image
  const colorScaleImage = document.createElementNS(NAME_SPACE_URI, "g");
  colorScaleImage.innerHTML = plasmaSvgString;
  colorScaleImage.setAttribute("transform", `translate(${xPosition}, 50)`);
  colorScaleImage.setAttribute("width", `${width}px`);

  // Color scale values
  const colorScaleValuesGroup = document.createElementNS(NAME_SPACE_URI, "g");
  colorScaleValuesGroup.setAttribute(
    "transform",
    `translate(${xPosition}, 80)`
  );
  const colorScaleValues = element.querySelectorAll(`.low-high span`);

  const colorScaleValueLow = document.createElementNS(NAME_SPACE_URI, "text");
  colorScaleValueLow.setAttribute("x", "0");
  colorScaleValueLow.textContent = colorScaleValues[0].innerHTML || null;

  const colorScaleValueHigh = document.createElementNS(NAME_SPACE_URI, "text");
  colorScaleValueHigh.setAttribute("x", `120`);
  colorScaleValueHigh.setAttribute("text-anchor", "end");
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
}: {
  heatmapContainer?: HTMLElement | null;
}) {
  if (!heatmapContainer) return;

  const xAxis = heatmapContainer.querySelector(`#x-axis-wrapper`);

  // Create root SVG elemnt
  const svg = document.createElementNS(NAME_SPACE_URI, "svg");

  const svgAttributes = {
    fill: ECHART_AXIS_LABEL_COLOR_HEX,
    "font-family": FONT_FAMILY,
    "font-size": ECHART_AXIS_LABEL_FONT_SIZE_PX,
    height: `${200}px`,
    x: `${Y_AXIS_CHART_WIDTH_PX}`,
    y: `10`,
  };

  applyAttributes(svg, svgAttributes);

  // Create cell type count label
  const cellTypeCount = document.createElementNS(NAME_SPACE_URI, "text");

  const cellTypeCountAttributes = {
    "text-anchor": "end",
    transform: `translate(35, ${X_AXIS_CHART_HEIGHT_PX}) rotate(90)`,
    width: `100px`,
    x: 0,
  };

  applyAttributes(cellTypeCount, cellTypeCountAttributes);
  cellTypeCount.textContent = "Cell Count";

  svg.append(cellTypeCount);

  // Create gene labels
  const geneLabelContainer = document.createElementNS(NAME_SPACE_URI, "g");

  Array.from(
    xAxis?.querySelectorAll(`[data-test-id*='gene-label-'] span`) || []
  ).forEach((label, index) => {
    const geneLabelText = document.createElementNS(NAME_SPACE_URI, "text");

    const labelAttributes = {
      "text-anchor": "end",
      transform: `translate(${
        65 + 20 * index
      }, ${X_AXIS_CHART_HEIGHT_PX}) rotate(90)`,
      x: 0,
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
  height,
}: {
  heatmapContainer?: HTMLElement | null;
  tissueName: Tissue;
  height: number;
}) {
  if (!heatmapContainer) return;

  const yAxis = heatmapContainer.querySelector(`#${tissueName}-y-axis`);

  const svg = document.createElementNS(NAME_SPACE_URI, "svg");

  const svgAttributes = {
    height: `${height}px`,
    width: `${Y_AXIS_CHART_WIDTH_PX + 60}px`, // Adds padding for current tissue label
    y: `${X_AXIS_CHART_HEIGHT_PX + 25}`,
  };

  applyAttributes(svg, svgAttributes);

  // Container for tissue label
  const tissueLabelContainer = document.createElementNS(NAME_SPACE_URI, "g");

  const tissueLabelText = document.createElementNS(NAME_SPACE_URI, "text");
  tissueLabelText.textContent = tissueName;
  tissueLabelText.setAttribute("transform", "translate(35, 100) rotate(270)");
  tissueLabelText.setAttribute("font-family", "Inter, sans-serif");
  tissueLabelText.setAttribute("font-size", "14px");
  tissueLabelText.setAttribute("font-weight", "bold");

  const tissueLabelRect = document.createElementNS(NAME_SPACE_URI, "rect");
  tissueLabelRect.setAttribute("x", "40");
  tissueLabelRect.setAttribute("width", "5");
  tissueLabelRect.setAttribute("height", height + "");

  tissueLabelContainer.append(tissueLabelText);
  tissueLabelContainer.append(tissueLabelRect);

  // cell type names container group
  const cellTypeNamesContainer = document.createElementNS(NAME_SPACE_URI, "g");

  Array.from(
    yAxis?.querySelectorAll("[data-test-id='cell-type-label-count']") || []
  ).forEach((labelCount, index) => {
    const label = labelCount.querySelector(
      "[data-test-id='cell-type-name']"
    )?.textContent;

    const count = labelCount.querySelector(
      "[data-test-id='cell-count']"
    )?.textContent;

    // group
    const group = document.createElementNS(NAME_SPACE_URI, "g");

    const groupAttributes = {
      fill: ECHART_AXIS_LABEL_COLOR_HEX,
      "font-family": FONT_FAMILY,
      "font-size": ECHART_AXIS_LABEL_FONT_SIZE_PX,
      /**
       * (thuang): Add `HEAT_MAP_BASE_CELL_PX / 2` top margin, so we render the
       * first label in the middle of the first cell
       */
      transform: `translate(50, ${
        HEAT_MAP_BASE_CELL_PX / 2 + index * HEAT_MAP_BASE_CELL_PX
      })`,
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
      y: 0,
    };

    applyAttributes(cellTypeLabel, labelAttributes);
    cellTypeLabel.textContent = String(label);

    // cellCount
    const cellCount = document.createElementNS(NAME_SPACE_URI, "text");

    const countAttributes = {
      "text-anchor": "end",
      x: Y_AXIS_CHART_WIDTH_PX,
      y: 0,
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

function applyAttributes(
  node: SVGElement,
  attributes: Record<string, string | number>
) {
  for (const [key, value] of Object.entries(attributes)) {
    node.setAttribute(key, String(value));
  }
}
