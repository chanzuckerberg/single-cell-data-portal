import { toPng, toSvg } from "html-to-image";
import { debounce } from "lodash";
import { Dispatch, useCallback, useContext, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { stringify as csvStringify } from "csv-stringify/sync";
import { HEATMAP_CONTAINER_ID } from "src/views/WheresMyGene/common/constants";
import { CellType, Tissue } from "src/views/WheresMyGene/common/types";
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
import {
  ButtonWrapper,
  DownloadButton,
  StyledButtonContainer,
  StyledDiv,
  StyledFormControlLabel,
  StyledInputCheckBox,
  StyledModal,
  StyledSection,
  StyledTitle,
} from "./style";
import {
  NAME_SPACE_URI,
  renderLegend,
  renderXAxis,
  renderYAxis,
} from "./utils";
import {
  FilterDimensions,
  OntologyTerm,
  SelectedFilters,
} from "src/common/queries/wheresMyGene";
import { StateContext } from "src/views/WheresMyGene/common/store";

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
  setDownloadStatus: React.Dispatch<
    React.SetStateAction<{
      isLoading: boolean;
      blur?: boolean;
    }>
  >;
  setEchartsRendererMode: Dispatch<React.SetStateAction<"canvas" | "svg">>;
  allChartProps: { [tissue: string]: ChartProps };
  availableFilters: Partial<FilterDimensions>;
  availableOrganisms: OntologyTerm[];
}

interface ExportData {
  input: string | ArrayBuffer;
  name: string;
}

export default function SaveImage({
  selectedTissues,
  selectedGenes,
  selectedCellTypes,
  setDownloadStatus,
  setEchartsRendererMode,
  allChartProps,
  availableFilters,
  availableOrganisms,
}: Props): JSX.Element {
  const { selectedFilters, selectedOrganismId } = useContext(StateContext);
  const [isOpen, setIsOpen] = useState(false);
  const [selectedFileTypes, setFileTypes] = useState<("png" | "svg" | "csv")[]>(
    ["png"]
  );
  const handleButtonClick = useCallback(() => {
    if (!isOpen) track(EVENTS.WMG_DOWNLOAD_CLICKED);
    setIsOpen(!isOpen);
  }, [isOpen]);

  function selectFileType(fileType: "png" | "svg" | "csv") {
    const index = selectedFileTypes.indexOf(fileType);
    if (index >= 0) {
      selectedFileTypes.splice(index, 1);
    } else {
      selectedFileTypes.push(fileType);
    }
    setFileTypes([...selectedFileTypes]);
  }

  const handleDownload = useCallback(async () => {
    setIsOpen(false);

    setDownloadStatus({ isLoading: true, blur: true });

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
        availableFilters,
        availableOrganisms,
        heatmapNode,
        observer,
        selectedCellTypes,
        selectedFileTypes,
        selectedFilters,
        selectedGenes,
        selectedOrganismId,
        selectedTissues,
        setEchartsRendererMode,
        setDownloadStatus,
      }),
      MUTATION_OBSERVER_TIMEOUT
    );

    observer.observe(heatmapNode, config);

    if (selectedFileTypes.includes("svg")) {
      setEchartsRendererMode("svg");
    }
  }, [
    allChartProps,
    availableFilters,
    availableOrganisms,
    selectedCellTypes,
    selectedFileTypes,
    selectedFilters,
    selectedGenes,
    selectedOrganismId,
    selectedTissues,
    setEchartsRendererMode,
    setDownloadStatus,
  ]);

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
      <StyledModal
        isOpen={isOpen}
        title="Download"
        onClose={handleButtonClick}
        isCloseButtonShown={false}
      >
        <div>
          <StyledSection>
            <StyledTitle>FIGURE</StyledTitle>
            <StyledDiv>
              <StyledFormControlLabel
                control={
                  <StyledInputCheckBox
                    checked={selectedFileTypes.includes("png")}
                  />
                }
                label="PNG"
                onChange={() => selectFileType("png")}
              />

              <StyledFormControlLabel
                control={
                  <StyledInputCheckBox
                    checked={selectedFileTypes.includes("svg")}
                  />
                }
                label="SVG"
                onChange={() => selectFileType("svg")}
              />
            </StyledDiv>
          </StyledSection>
          <StyledSection>
            <StyledTitle>DATA</StyledTitle>
            <StyledDiv>
              <StyledFormControlLabel
                control={
                  <StyledInputCheckBox
                    checked={selectedFileTypes.includes("csv")}
                  />
                }
                label="CSV"
                onChange={() => selectFileType("csv")}
              />
            </StyledDiv>
          </StyledSection>
        </div>

        <StyledButtonContainer>
          <DownloadButton
            sdsType="secondary"
            sdsStyle="minimal"
            onClick={handleButtonClick}
          >
            Cancel
          </DownloadButton>
          <DownloadButton
            sdsStyle="square"
            sdsType="primary"
            sdsSize="large"
            onClick={handleDownload}
            disabled={!selectedFileTypes.length}
          >
            Download
          </DownloadButton>
        </StyledButtonContainer>
      </StyledModal>
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
  selectedGenes: Props["selectedGenes"],
  tissue: string,
  availableFilters: Partial<FilterDimensions>,
  selectedFilters: SelectedFilters,
  selectedOrganismId: string | null,
  availableOrganisms: OntologyTerm[]
) {
  const { datasets, disease_terms, self_reported_ethnicity_terms, sex_terms } =
    availableFilters;

  const output: (string | number)[][] = [];

  // Metadata as comments
  output.push([`# ${new Date().toString()}`]);

  // Dataset
  output.push(["# Dataset"]);
  output.push([
    `${
      datasets
        ?.filter((option) => {
          return selectedFilters.datasets.includes(option.id);
        })
        .map((selected) => selected.label)
        .join(", ") || ""
    }`,
  ]);

  // Disease
  output.push(["# Disease"]);
  output.push([
    `${
      disease_terms
        ?.filter((option) => {
          return selectedFilters.diseases.includes(option.id);
        })
        .map((selected) => selected.name)
        .join(", ") || ""
    }`,
  ]);

  // Organism
  output.push(["# Self-Reported Ethnicity"]);
  output.push([
    `${
      self_reported_ethnicity_terms
        ?.filter((option) => {
          return selectedFilters.ethnicities.includes(option.id);
        })
        .map((selected) => selected.name)
        .join(", ") || ""
    }`,
  ]);

  // Sex
  output.push(["# Sex"]);
  output.push([
    `${
      sex_terms
        ?.filter((option) => {
          return selectedFilters.sexes.includes(option.id);
        })
        .map((selected) => selected.name)
        .join(", ") || ""
    }`,
  ]);

  // Organism
  output.push(["# Organism"]);
  output.push([
    `${
      availableOrganisms.find((organism) => organism.id === selectedOrganismId)
        ?.name
    }`,
  ]);

  // Column Names
  output.push([
    "Tissue",
    "Cell Type",
    "Cell Count",
    "Tissue Composition",
    "Gene Symbol",
    "Expression",
    "Expression, Scaled",
    "Number of Cells Expressing Genes",
  ]);

  for (const cellTypeMetaData of allChartProps[
    tissue
  ].cellTypeMetadata.reverse()) {
    const { id, name, total_count } =
      deserializeCellTypeMetadata(cellTypeMetaData);

    for (const geneName of selectedGenes) {
      const geneExpression = allChartProps[tissue].chartData.find(
        (value) => value.id === `${id}-${geneName}`
      );

      output.push([
        tissue,
        name,
        total_count,
        Number((geneExpression?.tissuePercentage || 0) * 100).toFixed(2) + "%",
        geneName,
        geneExpression?.meanExpression ?? "",
        geneExpression?.scaledMeanExpression ?? "",
        geneExpression?.expressedCellCount ?? "",
      ]);
    }
  }

  return csvStringify(output);
}

async function generateImage(
  fileType: string,
  heatmapNode: HTMLDivElement,
  height: number,
  tissueName: string,
  isMultipleTissues = false
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
async function initiateDownload(
  selectedFileTypes: string[],
  exports: ExportData[]
) {
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
    } else {
      const encodedContent = encodeURIComponent(exports[0].input as string);

      if (fileType === "svg") {
        // (ashin-czi): Fixes SVG string breaking when encountering a "#" character
        link.href = `data:image/svg+xml,${encodedContent}`;
      } else if (fileType === "csv") {
        link.href = `data:text/csv,${encodedContent}`;
      }
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
  availableFilters,
  availableOrganisms,
  selectedOrganismId,
  selectedFilters,
  observer,
  setDownloadStatus,
  setEchartsRendererMode,
}: {
  allChartProps: { [tissue: string]: ChartProps };
  heatmapNode: HTMLDivElement;
  selectedTissues: Props["selectedTissues"];
  selectedGenes: Props["selectedGenes"];
  selectedCellTypes: Props["selectedCellTypes"];
  selectedFileTypes: ("png" | "svg" | "csv")[];
  availableFilters: Partial<FilterDimensions>;
  availableOrganisms: OntologyTerm[];
  selectedOrganismId: string | null;
  selectedFilters: SelectedFilters;
  observer: MutationObserver;
  setDownloadStatus: React.Dispatch<
    React.SetStateAction<{
      isLoading: boolean;
      blur?: boolean;
    }>
  >;
  setEchartsRendererMode: (mode: "canvas" | "svg") => void;
}) {
  return async () => {
    try {
      const isPng = selectedFileTypes.includes("png");

      //(ashin): #3569 Get scrollTop to go back to place after downloading image
      const heatmapContainer = document.getElementById(
        HEATMAP_CONTAINER_ID
      ) as HTMLCanvasElement;
      heatmapContainerScrollTop = heatmapContainer?.scrollTop;

      const initialWidth = heatmapNode.style.width;

      if (isPng) {
        // Adding this class causes the y-axis scrolling to jump but is required for image download
        heatmapNode.classList.add("CLONED");

        heatmapNode.style.width = `${
          getHeatmapWidth(selectedGenes) + Y_AXIS_CHART_WIDTH_PX + 100
        }px`;
      }

      const exports: ExportData[] = (
        await Promise.all(
          selectedTissues.map(async (tissueName) => {
            // Handles if whitespace is in the tissue name for the element ID
            const formattedTissueName = tissueName.replace(/\s+/g, "-");

            // Generate exports for each filetype for the tissue
            return await Promise.all(
              selectedFileTypes.map(async (fileType) => {
                let input;

                if (fileType === "csv") {
                  input = generateCsv(
                    allChartProps,
                    selectedGenes,
                    tissueName,
                    availableFilters,
                    selectedFilters,
                    selectedOrganismId,
                    availableOrganisms
                  );
                } else {
                  const height = getHeatmapHeight(
                    selectedCellTypes[tissueName]
                  );

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
              })
            );
          })
        )
      ).flat();

      if (isPng) {
        //(thuang): #3569 Restore scrollTop position
        heatmapNode.classList.remove("CLONED");
        heatmapNode.style.width = initialWidth;
        if (heatmapContainer) {
          heatmapContainer.scrollTop = heatmapContainerScrollTop || 0;
        }
      }

      await initiateDownload(selectedFileTypes, exports);

      track(EVENTS.WMG_DOWNLOAD_COMPLETE, {
        file_type: selectedFileTypes.join(","),

        dataset_filter: selectedFilters.datasets.join(","),
        disease_filter: selectedFilters.diseases.join(","),
        self_reported_ethnicity_filter: selectedFilters.ethnicities.join(","),
        sex_filter: selectedFilters.sexes.join(","),
        group_by_option: "", // TODO: AFTER COMPARE

        genes: selectedGenes.toString(),
        tissues: selectedTissues.toString(),
      });
    } catch (error) {
      console.error(error);
    }

    observer.disconnect();
    setEchartsRendererMode("canvas");
    setDownloadStatus({ isLoading: false });
  };
}
