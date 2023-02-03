import { Classes, Intent } from "@blueprintjs/core";
import { FormControlLabel } from "@mui/material";
import { InputRadio } from "czifui";
import { toPng, toSvg } from "html-to-image";
import { Dispatch, useCallback, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  Section,
  Title,
} from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/common/style";
import Modal from "src/components/common/Modal";
import { HEATMAP_CONTAINER_ID } from "src/views/WheresMyGene/common/constants";
import { CellType, Tissue } from "src/views/WheresMyGene/common/types";
import {
  ECHART_AXIS_LABEL_COLOR_HEX,
  ECHART_AXIS_LABEL_FONT_SIZE_PX,
} from "../../../HeatMap/components/XAxisChart/style";
import {
  getHeatmapHeight,
  getHeatmapWidth,
  HEAT_MAP_BASE_CELL_PX,
  X_AXIS_CHART_HEIGHT_PX,
  Y_AXIS_CHART_WIDTH_PX,
} from "../../../HeatMap/utils";
import { plasmaSvgString } from "../../../InfoPanel/components/ColorScale";
import { Label } from "../../style";
import { StyledButtonIcon } from "../QuickSelect/style";
import { ButtonWrapper, DownloadButton, StyledDiv } from "./style";

let heatmapContainerScrollTop: number | undefined;

export const EXCLUDE_IN_SCREENSHOT_CLASS_NAME = "screenshot-exclude";
const screenshotFilter =
  (tissue: string) =>
  (domNode: HTMLElement): boolean => {
    const isScreenshotJank = domNode.classList?.contains(
      EXCLUDE_IN_SCREENSHOT_CLASS_NAME
    );
    const isNoScript = domNode.tagName === "NOSCRIPT";

    const isNonTissueChart =
      domNode.id &&
      domNode.id.includes("chart") &&
      !domNode.id.includes(tissue);
    const isNonTissueYAxis =
      domNode.id &&
      domNode.id.includes("y-axis") &&
      !domNode.id.includes(tissue);
    return (
      !isScreenshotJank && !isNoScript && !isNonTissueChart && !isNonTissueYAxis
    );
  };

function base64URLToArrayBuffer(url: string) {
  // parse dataURL to base64 string
  const base64 = url.split(",")[1];
  const binary_string = window.atob(base64);
  const len = binary_string.length;
  const bytes = new Uint8Array(len);
  for (let i = 0; i < len; i++) {
    bytes[i] = binary_string.charCodeAt(i);
  }
  return bytes.buffer;
}

interface Props {
  selectedTissues: Array<string>;
  selectedGenes: Array<string>;
  selectedCellTypes: { [tissue: string]: CellType[] };
  setIsDownloading: (isDownloading: boolean) => void;
  setEchartsRendererMode: Dispatch<React.SetStateAction<"canvas" | "svg">>;
}

export default function SaveImage({
  selectedTissues,
  selectedGenes,
  selectedCellTypes,
  setIsDownloading,
  setEchartsRendererMode,
}: Props): JSX.Element {
  const [isOpen, setIsOpen] = useState(false);
  const [fileType, setFileType] = useState<"png" | "svg">("png");
  const handleButtonClick = useCallback(() => {
    if (!isOpen) track(EVENTS.WMG_DOWNLOAD_CLICKED);
    setIsOpen(!isOpen);
  }, [isOpen]);

  const handleDownload = useCallback(async () => {
    setIsOpen(false);
    setIsDownloading(true);
    setEchartsRendererMode("svg");

    await new Promise((resolve) => setTimeout(resolve, 5 * 1000));

    try {
      const heatmapNode = document.getElementById("view") as HTMLCanvasElement;
      //(ashin): #3569 Get scrollTop to go back to place after downloading image
      const heatmapContainer = document.getElementById(
        HEATMAP_CONTAINER_ID
      ) as HTMLCanvasElement;
      heatmapContainerScrollTop = heatmapContainer?.scrollTop;

      // Adding this class causes the y-axis scrolling to jump but is required for image download
      heatmapNode.classList.add("CLONED");
      const initialWidth = heatmapNode.style.width;
      heatmapNode.style.width = `${
        getHeatmapWidth(selectedGenes) + Y_AXIS_CHART_WIDTH_PX + 100
      }px`;
      const isPNG = fileType === "png";
      const convertHTMLtoImage = isPNG ? toPng : toSvg;
      const images = await Promise.all(
        selectedTissues.map(async (tissueName) => {
          const height = getHeatmapHeight(selectedCellTypes[tissueName]);

          // Handles if whitespace is in the tissue name for the element ID
          tissueName = tissueName.replace(/\s+/g, "-");

          const imageURL = await convertHTMLtoImage(heatmapNode, {
            backgroundColor: "white",
            filter: screenshotFilter(tissueName),
            height: height + X_AXIS_CHART_HEIGHT_PX + 120,
            pixelRatio: 4,
            width: heatmapNode.width,
          });
          // raw URI if only one tissue is selected

          const input = isPNG
            ? selectedTissues.length === 1
              ? imageURL
              : base64URLToArrayBuffer(imageURL)
            : processSvg({
                height,
                svg: decodeURIComponent(imageURL.split(",")[1]),
                tissueName,
              });

          return {
            input,
            name: `${tissueName}.${fileType}`,
          };
        })
      );

      //(thuang): #3569 Restore scrollTop position
      heatmapNode.classList.remove("CLONED");
      heatmapNode.style.width = initialWidth;
      if (heatmapContainer) {
        heatmapContainer.scrollTop = heatmapContainerScrollTop || 0;
      }

      const link = document.createElement("a");

      if (images.length > 1) {
        const { downloadZip } = await import("client-zip");
        const blob = await downloadZip(images).blob();
        link.href = URL.createObjectURL(blob);
        link.download = `CELLxGENE_gene_expression.zip`;
      } else {
        if (isPNG) {
          // PNG download link
          link.href = images[0].input as string;
          link.download = images[0].name;
        } else {
          // SVG download link
          // (ashin-czi): Fixes SVG string breaking when encountering a "#" character
          link.href =
            "data:image/svg+xml," +
            (images[0].input as string).replace(/#/g, "%23");
          link.download = images[0].name;
        }
      }
      link.click();
      link.remove();
      track(EVENTS.WMG_DOWNLOAD_COMPLETE, {
        file_type: fileType,
        genes: selectedGenes.toString(),
        tissues: selectedTissues.toString(),
      });
    } catch (error) {
      console.error(error);
    }

    setIsDownloading(false);
    setEchartsRendererMode("canvas");
  }, [
    fileType,
    selectedCellTypes,
    selectedTissues,
    selectedGenes,
    setIsDownloading,
    setEchartsRendererMode,
  ]);

  return (
    <>
      <ButtonWrapper className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
        <Label>Download</Label>
        <StyledButtonIcon
          data-test-id="download-button"
          // TODO: put handleButtonClick when svgs are fixed
          onClick={handleButtonClick}
          disabled={selectedTissues.length === 0 || selectedGenes.length === 0}
          sdsSize="medium"
          sdsType="primary"
          sdsIcon="download"
        />
      </ButtonWrapper>
      <Modal
        isOpen={isOpen}
        title="Download Figure"
        onClose={handleButtonClick}
        isCloseButtonShown={false}
      >
        <div className={Classes.DIALOG_BODY}>
          <Section>
            <Title>Image Format</Title>
            <StyledDiv>
              <FormControlLabel
                control={
                  <InputRadio
                    stage={fileType === "png" ? "checked" : "unchecked"}
                  />
                }
                label="PNG"
                onChange={() => setFileType("png")}
              />

              <FormControlLabel
                control={
                  <InputRadio
                    stage={fileType === "svg" ? "checked" : "unchecked"}
                  />
                }
                label="SVG"
                onChange={() => setFileType("svg")}
              />
            </StyledDiv>
          </Section>
        </div>

        <div className={Classes.DIALOG_FOOTER}>
          <DownloadButton onClick={handleButtonClick} minimal>
            Cancel
          </DownloadButton>
          <DownloadButton intent={Intent.PRIMARY} onClick={handleDownload}>
            Download
          </DownloadButton>
        </div>
      </Modal>
    </>
  );
}

function processSvg({
  svg,
  tissueName,
  height,
}: {
  svg: string;
  tissueName: Tissue;
  height: number;
}) {
  const heatmapNode = new DOMParser().parseFromString(svg, "image/svg+xml");
  const heatmapContainer = heatmapNode
    .querySelector("foreignObject")
    ?.querySelector("div");

  // Render elements to SVG
  const xAxisSvg = renderXAxis({ heatmapContainer });
  const yAxisSvg = renderYAxis({ heatmapContainer, height, tissueName });
  const dotsSvg = renderDots({
    heatmapContainer,
    tissueName,
    xPosition: yAxisSvg!.width.baseVal.value,
  });
  const legendSvg = renderLegend({ heatmapContainer });

  const svgWidth =
    yAxisSvg!.width.baseVal.value + dotsSvg!.width.baseVal.value + 20;
  const svgHeight =
    dotsSvg!.height.baseVal.value + xAxisSvg!.height.baseVal.value;
  const finalSvg = document.createElementNS(NAME_SPACE_URI, "svg");

  finalSvg.setAttributeNS(
    "http://www.w3.org/2000/xmlns/",
    "xmlns",
    NAME_SPACE_URI
  );
  finalSvg.setAttribute("height", `${svgHeight}`);
  finalSvg.setAttribute("width", `${svgWidth}`);

  // Append elements to final SVG
  finalSvg.append(legendSvg!);
  finalSvg.append(xAxisSvg!);
  finalSvg.append(yAxisSvg!);
  finalSvg.append(dotsSvg!);

  return finalSvg.outerHTML;
}

function renderDots({
  heatmapContainer,
  tissueName,
  xPosition,
}: {
  heatmapContainer?: HTMLElement | null;
  tissueName: Tissue;
  xPosition: number;
}) {
  if (!heatmapContainer) return;

  const chart = heatmapContainer
    .querySelector(`#${tissueName}-chart`)
    ?.querySelector("svg");

  chart?.setAttribute("y", `${X_AXIS_CHART_HEIGHT_PX + 20}`);
  chart?.setAttribute("x", `${xPosition}`);

  // Cleanup as style attributes aren't used in SVG files
  chart?.removeAttribute("style");
  Array.from(chart?.querySelectorAll("rect, path, g") || []).forEach(
    (element) => {
      element.removeAttribute("style");
    }
  );

  return chart;
}

const NAME_SPACE_URI = "http://www.w3.org/2000/svg";

function renderLegend({
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

  const FONT_FAMILY = "sans-serif";

  // Create root SVG elemnt
  const svg = document.createElementNS(NAME_SPACE_URI, "svg");

  const svgAttributes = {
    x: `0`,
    y: `0`,
    fill: ECHART_AXIS_LABEL_COLOR_HEX,
    "font-family": FONT_FAMILY,
    "font-size": ECHART_AXIS_LABEL_FONT_SIZE_PX,
  };

  applyAttributes(svg, svgAttributes);

  svg.append(colorScale!);
  svg.append(expressedInCells!);

  return svg;
}

function renderExpressedInCells({
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
  expressedCellsLabel.textContent = element.querySelector(
    `#expressed-in-cells-label`
  )?.innerHTML!;
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
  expressedCellsValueLow.textContent = expressedCellsValues[0].innerHTML!;

  const expressedCellsValueHigh = document.createElementNS(
    NAME_SPACE_URI,
    "text"
  );
  expressedCellsValueHigh.setAttribute("x", `${width}`);
  expressedCellsValueHigh.setAttribute("text-anchor", "end");
  expressedCellsValueHigh.textContent = expressedCellsValues[1].innerHTML!;

  expressedCellsValuesGroup.append(expressedCellsValueLow);
  expressedCellsValuesGroup.append(expressedCellsValueHigh);

  // Append color scale legend
  expressedCellsGroup.append(expressedCellsLabel);
  expressedCellsGroup.append(expressedCellsDots);
  expressedCellsGroup.append(expressedCellsValuesGroup);

  return expressedCellsGroup;
}

function renderColorScale({
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
  colorScaleLabel.textContent = element.querySelector(
    `#relative-gene-expression-label`
  )?.innerHTML!;
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
  colorScaleValueLow.textContent = colorScaleValues[0].innerHTML!;

  const colorScaleValueHigh = document.createElementNS(NAME_SPACE_URI, "text");
  colorScaleValueHigh.setAttribute("x", `120`);
  colorScaleValueHigh.setAttribute("text-anchor", "end");
  colorScaleValueHigh.textContent = colorScaleValues[1].innerHTML!;

  colorScaleValuesGroup.append(colorScaleValueLow);
  colorScaleValuesGroup.append(colorScaleValueHigh);

  // Append color scale legend
  colorScaleGroup.append(colorScaleLabel);
  colorScaleGroup.append(colorScaleImage);
  colorScaleGroup.append(colorScaleValuesGroup);

  return colorScaleGroup;
}

function renderXAxis({
  heatmapContainer,
}: {
  heatmapContainer?: HTMLElement | null;
}) {
  if (!heatmapContainer) return;

  const xAxis = heatmapContainer.querySelector(`#x-axis-wrapper`);

  const FONT_FAMILY = "sans-serif";

  // Create root SVG elemnt
  const svg = document.createElementNS(NAME_SPACE_URI, "svg");

  const svgAttributes = {
    height: `${200}px`,
    x: `${Y_AXIS_CHART_WIDTH_PX}`,
    y: `10`,
    fill: ECHART_AXIS_LABEL_COLOR_HEX,
    "font-family": FONT_FAMILY,
    "font-size": ECHART_AXIS_LABEL_FONT_SIZE_PX,
  };

  applyAttributes(svg, svgAttributes);

  // Create cell type count label
  const cellTypeCount = document.createElementNS(NAME_SPACE_URI, "text");

  const cellTypeCountAttributes = {
    x: 0,
    width: `100px`,
    transform: `translate(40, ${X_AXIS_CHART_HEIGHT_PX}) rotate(90)`,
    "text-anchor": "end",
  };

  applyAttributes(cellTypeCount, cellTypeCountAttributes);
  cellTypeCount.textContent = "Cell Count";

  svg.append(cellTypeCount);

  // Create gene labels
  const geneLabelContainer = document.createElementNS(NAME_SPACE_URI, "g");

  Array.from(
    xAxis?.querySelectorAll(`[data-test-id*='gene-label-'] div`) || []
  ).forEach((label, index) => {
    const geneLabelText = document.createElementNS(NAME_SPACE_URI, "text");

    const labelAttributes = {
      x: 0,
      transform: `translate(${
        65 + 20 * index
      }, ${X_AXIS_CHART_HEIGHT_PX}) rotate(90)`,
      "text-anchor": "end",
    };

    applyAttributes(geneLabelText, labelAttributes);
    geneLabelText.textContent = label.innerHTML;

    geneLabelContainer.append(geneLabelText);
  });

  svg.append(geneLabelContainer);

  return svg;
}

function renderYAxis({
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

  const FONT_FAMILY = "sans-serif";

  const svg = document.createElementNS(NAME_SPACE_URI, "svg");

  const svgAttributes = {
    y: `${X_AXIS_CHART_HEIGHT_PX + 25}`,
    height: `${height}px`,
    width: `${Y_AXIS_CHART_WIDTH_PX + 60}px`, // Adds padding for current tissue label
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
      "[data-test-id='cell-type-label']"
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
