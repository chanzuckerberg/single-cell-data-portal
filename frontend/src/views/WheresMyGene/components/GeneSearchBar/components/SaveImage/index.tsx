/* eslint-disable prettier/prettier */
import { Classes, Intent } from "@blueprintjs/core";
import { FormControlLabel } from "@mui/material";
import { InputCheckbox } from "czifui";
import { toPng, toSvg } from "html-to-image";
import { debounce } from "lodash";
import { Dispatch, useCallback, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  Section,
  Title,
} from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/common/style";
import Modal from "src/components/common/Modal";
import { HEATMAP_CONTAINER_ID } from "src/views/WheresMyGene/common/constants";
import {
  CellType,
  Tissue,
} from "src/views/WheresMyGene/common/types";
import { ChartProps } from "../../../HeatMap/hooks/common/types";

import {
  deserializeCellTypeMetadata,
  getHeatmapHeight,
  getHeatmapWidth,
  X_AXIS_CHART_HEIGHT_PX,
  Y_AXIS_CHART_WIDTH_PX,
} from "../../../HeatMap/utils";

import { Label } from "../../style";
import { StyledButtonIcon } from "../QuickSelect/style";
import { ButtonWrapper, DownloadButton, StyledDiv } from "./style";
import {
  NAME_SPACE_URI,
  renderLegend,
  renderXAxis,
  renderYAxis,
} from "./utils";

let heatmapContainerScrollTop: number | undefined;

const MUTATION_OBSERVER_TIMEOUT = 3 * 1000;

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
  allChartProps: { [tissue: string]: ChartProps };
}

interface ExportData {
  input: string | ArrayBuffer;
  name: string;
}

export default function SaveImage({
  selectedTissues,
  selectedGenes,
  selectedCellTypes,
  setIsDownloading,
  setEchartsRendererMode,
  allChartProps,
}: Props): JSX.Element {
  const [isOpen, setIsOpen] = useState(false);
  const [selectedFileTypes, setFileTypes] = useState<("png" | "svg" | "csv")[]>(["png"]);
  const handleButtonClick = useCallback(() => {
    if (!isOpen) track(EVENTS.WMG_DOWNLOAD_CLICKED);
    setIsOpen(!isOpen);
  }, [isOpen]);

  function selectFileType(fileType: "png" | "svg" | "csv") {
    const index = selectedFileTypes.indexOf(fileType);
    if(index >= 0) {
      selectedFileTypes.splice(index, 1);
    } else {
      selectedFileTypes.push(fileType);
    }
    console.log(selectedFileTypes);
    setFileTypes([...selectedFileTypes]);
  }

  const handleDownload = useCallback(async () => {
    setIsOpen(false);
    setIsDownloading(true);

    const heatmapNode = document.getElementById("view") as HTMLDivElement;

    // Options for the observer (which mutations to observe)
    const config = { childList: true, subtree: true };

    // Callback function to execute when mutations are observed
    const callback = debounce(() => {
      download();
    }, MUTATION_OBSERVER_TIMEOUT);

    const observer = new MutationObserver(callback);

    /**
     * (thuang): We only want to call download() once, AFTER the heatmap is fully rendered.
     * We assume the heatmap is fully rendered after no more changes to the DOM tree
     * after `MUTATION_OBSERVER_TIMEOUT` milliseconds
     */
    const download = debounce(
      download_({
        allChartProps,
        heatmapNode,
        observer,
        selectedCellTypes,
        selectedFileTypes,
        selectedGenes,
        selectedTissues,
        setEchartsRendererMode,
        setIsDownloading,
      }),
      MUTATION_OBSERVER_TIMEOUT
    );

    observer.observe(heatmapNode, config);

    setEchartsRendererMode("svg");
  }, [setIsDownloading, allChartProps, selectedCellTypes, selectedFileTypes, selectedGenes, selectedTissues, setEchartsRendererMode]);

  return (
    <>
      <ButtonWrapper className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
        <Label>Download</Label>
        <StyledButtonIcon
          data-test-id="download-button"
          onClick={handleButtonClick}
          disabled={selectedTissues.length === 0 || selectedGenes.length === 0}
          sdsSize="medium"
          sdsType="primary"
          sdsIcon="download"
        />
      </ButtonWrapper>
      <Modal
        isOpen={isOpen}
        title="DOWNLOAD"
        onClose={handleButtonClick}
        isCloseButtonShown={false}
      >
        <div className={Classes.DIALOG_BODY}>
          <Section>
            <Title>FIGURE</Title>
            <StyledDiv>
              <FormControlLabel
                control={
                  <InputCheckbox 
                    checked={selectedFileTypes.includes("png")}
                  />
                }
                label="PNG"
                onChange={() => selectFileType("png")}
              />

              <FormControlLabel
                control={
                  <InputCheckbox 
                    checked={selectedFileTypes.includes("svg")}
                  />
                }
                label="SVG"
                onChange={() => selectFileType("svg")}
              />
            </StyledDiv>
          </Section>
          <Section>
            <Title>DATA</Title>
            <StyledDiv>
              <FormControlLabel
                control={
                  <InputCheckbox 
                    checked={selectedFileTypes.includes("csv")}
                  />
                }
                label="CSV"
                onChange={() => selectFileType("csv")}
              />
            </StyledDiv>
          </Section>
        </div>

        <div className={Classes.DIALOG_FOOTER}>
          <DownloadButton onClick={handleButtonClick} minimal>
            Cancel
          </DownloadButton>
          <DownloadButton intent={Intent.PRIMARY} onClick={handleDownload} disabled={!selectedFileTypes.length}>
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
  finalSvg.append(legendSvg || "");
  finalSvg.append(xAxisSvg || "");
  finalSvg.append(yAxisSvg || "");
  finalSvg.append(dotsSvg || "");

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

/**
 * Generates a CSV string for the tissue, cell type, and gene expression combinations
 * @param selectedCellTypes 
 * @param selectedGenes 
 * @returns {string}
 */
function generateCsv(
  allChartProps: { [tissue: string]: ChartProps },
  selectedGenes:  Props["selectedGenes"],
  tissue: string,
) {
  
  let rows: string = 
  [
    '"Tissue"',
    '"Cell Type"',
    '"Cell Count"',
    '"Tissue Composition"',
    '"Gene Symbol"',
    '"Expression"',
    '"Expression, Scaled"',
    '"Number of Cells Expressing Genes"\n',
  ].join();
  
  console.log(JSON.stringify(allChartProps, null, 4));

  for (const cellTypeMetaData of allChartProps[tissue].cellTypeMetadata.reverse()) {
    const { 
      id, 
      name, 
      tissue, 
      total_count 
    } = deserializeCellTypeMetadata(cellTypeMetaData);

    for (const geneName of selectedGenes) {
      const geneExpression = allChartProps[tissue].chartData.find((value) => value.id === `${id}-${geneName}`);

      // console.log("expression " + JSON.stringify(geneExpression))

      rows += (`"${[
        tissue,
        name,
        total_count,
        Number((geneExpression?.tissuePercentage || 0) * 100).toFixed(2) + "%",
        geneName,
        geneExpression?.meanExpression ?? "",
        geneExpression?.scaledMeanExpression ?? "",
        geneExpression?.expressedCellCount ?? "",
      ].join('","')}"\n`);
    }
  }

  return rows;
}

async function generateImage(
  fileType: string,
  heatmapNode: HTMLDivElement,
  height: number,
  tissueName: string,
  isMultipleTissues = false,
): Promise<string | ArrayBuffer> {
  const convertHTMLtoImage = fileType === "png" ? toPng : toSvg;
  
  const imageURL = await convertHTMLtoImage(heatmapNode, {
    backgroundColor: "white",
    filter: screenshotFilter(tissueName),
    height: height + X_AXIS_CHART_HEIGHT_PX + 120,
    pixelRatio: 4,
    width: heatmapNode.offsetWidth,
  });

  let input: string | ArrayBuffer = imageURL;

  if (fileType === "svg") {
    input = processSvg({
      height,
      svg: decodeURIComponent(imageURL.split(",")[1]),
      tissueName,
    });
  } else if (fileType === "png" && isMultipleTissues) {
    input = base64URLToArrayBuffer(imageURL);
  }

  return input;
}

/**
 * Generates an anchor element, populates the href and download properties, and clicks the element to download. 
 * Element is deleted after clicked
 * @param selectedFileTypes 
 * @param exports 
 */
async function initiateDownload(selectedFileTypes: string[], exports: ExportData[]) {
  const link = document.createElement("a");

  if (exports.length > 1) {
    const { downloadZip } = await import("client-zip");
    const blob = await downloadZip(exports).blob();

    link.download = `CELLxGENE_gene_expression.zip`;
    link.href = URL.createObjectURL(blob);
  } else {
    // If there's only one export then there will be one filetype selected
    const fileType = selectedFileTypes[0];

    link.download = exports[0].name;

    if (fileType === "png") {
      link.href = exports[0].input as string;
    } else if (fileType === "svg") {
      // (ashin-czi): Fixes SVG string breaking when encountering a "#" character
      link.href = `data:image/svg+xml,${(exports[0].input as string).replace(/#/g, "%23")}`;
    } else if (fileType === "csv") {
      link.href = `data:text/csv,${exports[0].input as string}`;
    }
  }

  link.click();
  link.remove();
}

function download_({
  allChartProps,
  heatmapNode,
  selectedTissues,
  selectedGenes,
  selectedCellTypes,
  selectedFileTypes,
  observer,
  setIsDownloading,
  setEchartsRendererMode,
}: {
  allChartProps: { [tissue: string]: ChartProps };
  heatmapNode: HTMLDivElement;
  selectedTissues: Props["selectedTissues"];
  selectedGenes: Props["selectedGenes"];
  selectedCellTypes: Props["selectedCellTypes"];
  selectedFileTypes: ("png" | "svg" | "csv")[];
  observer: MutationObserver;
  setIsDownloading: (isDownloading: boolean) => void;
  setEchartsRendererMode: (mode: "canvas" | "svg") => void;
}) {
  return async () => {
    try {
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

      const exports: ExportData[] = (await Promise.all(
        selectedTissues.map(async (tissueName) => {
          // Handles if whitespace is in the tissue name for the element ID
          const formattedTissueName = tissueName.replace(/\s+/g, "-");

          // Generate exports for each filetype for the tissue
          return await Promise.all(selectedFileTypes.map(async (fileType) => {
            let input; 

            if (fileType === "csv") {
              input = generateCsv(allChartProps, selectedGenes, tissueName);
              // console.log(input)
            } else {
              const height = getHeatmapHeight(selectedCellTypes[tissueName]);

              input = await generateImage(
                fileType,
                heatmapNode,
                height,
                formattedTissueName,
                selectedTissues.length > 1 || selectedFileTypes.length > 1
              );
            }
  
            return {
              input,
              name: `${formattedTissueName}.${fileType}`,
            };
          }))
        })
      )).flat();

      
      console.log(exports);

      //(thuang): #3569 Restore scrollTop position
      heatmapNode.classList.remove("CLONED");
      heatmapNode.style.width = initialWidth;
      if (heatmapContainer) {
        heatmapContainer.scrollTop = heatmapContainerScrollTop || 0;
      }
      
      await initiateDownload(selectedFileTypes, exports);

      track(EVENTS.WMG_DOWNLOAD_COMPLETE, {
        file_type: selectedFileTypes.join(","),
        genes: selectedGenes.toString(),
        tissues: selectedTissues.toString(),
      });
    } catch (error) {
      console.error(error);
    }

    observer.disconnect();
    setIsDownloading(false);
    setEchartsRendererMode("canvas");
  };
}
