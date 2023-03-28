import { toPng, toSvg } from "html-to-image";
import { debounce } from "lodash";
import {
  Dispatch,
  SetStateAction,
  useCallback,
  useContext,
  useState,
} from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { stringify as csvStringify } from "csv-stringify/sync";
import {
  CompareId,
  getCompareOptionNameById,
  HEATMAP_CONTAINER_ID,
} from "src/views/WheresMyGene/common/constants";
import { CellType, Tissue } from "src/views/WheresMyGene/common/types";
import { ChartProps } from "../../../HeatMap/hooks/common/types";

import {
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
  StyledModalContent,
  StyledButtonContainer,
  StyledInputCheckboxWrapper,
  StyledModal,
  StyledSection,
  StyledTitle,
  StyledMessage,
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
  useAvailableOrganisms,
} from "src/common/queries/wheresMyGene";
import { State, StateContext } from "src/views/WheresMyGene/common/store";
import {
  buildCellTypeIdToMetadataMapping,
  csvGeneExpressionRow,
  csvHeaders,
} from "./csvUtils";
import { InputCheckbox } from "czifui";

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

export interface Props {
  selectedTissues: Array<string>;
  selectedGenes: Array<string>;
  selectedCellTypes: { [tissue: string]: CellType[] };
  setDownloadStatus: Dispatch<
    SetStateAction<{
      isLoading: boolean;
      blur?: boolean;
    }>
  >;
  setEchartsRendererMode: Dispatch<SetStateAction<"canvas" | "svg">>;
  allChartProps: { [tissue: string]: ChartProps };
  availableFilters: Partial<FilterDimensions>;
}

interface ExportData {
  input: string | ArrayBuffer;
  name: string;
}

export default function SaveExport({
  selectedTissues,
  selectedGenes,
  selectedCellTypes,
  setDownloadStatus,
  setEchartsRendererMode,
  allChartProps,
  availableFilters,
}: Props): JSX.Element {
  const { selectedFilters, selectedOrganismId, compare } =
    useContext(StateContext);
  const [isOpen, setIsOpen] = useState(false);
  const [selectedFileTypes, setFileTypes] = useState<("png" | "svg" | "csv")[]>(
    ["png"]
  );
  const { data: availableOrganisms } = useAvailableOrganisms();

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
        compare,
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
    setDownloadStatus,
    allChartProps,
    availableFilters,
    availableOrganisms,
    compare,
    selectedCellTypes,
    selectedFileTypes,
    selectedFilters,
    selectedGenes,
    selectedOrganismId,
    selectedTissues,
    setEchartsRendererMode,
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
        <StyledModalContent>
          <StyledSection>
            <StyledTitle>FIGURE</StyledTitle>
            <StyledInputCheckboxWrapper>
              <InputCheckbox
                label="PNG"
                caption="PNG CAPTION"
                onChange={() => selectFileType("png")}
                data-test-id="png-checkbox"
                checked={selectedFileTypes.includes("png")}
              />
            </StyledInputCheckboxWrapper>
            <StyledInputCheckboxWrapper>
              <InputCheckbox
                label="SVG"
                caption="SVG CAPTION"
                onChange={() => selectFileType("svg")}
                data-test-id="svg-checkbox"
                checked={selectedFileTypes.includes("svg")}
              />
            </StyledInputCheckboxWrapper>
          </StyledSection>

          <StyledSection>
            <StyledTitle>DATA</StyledTitle>
            <StyledInputCheckboxWrapper>
              <InputCheckbox
                label="CSV"
                caption="CSV CAPTION"
                onChange={() => selectFileType("csv")}
                data-test-id="csv-checkbox"
                checked={selectedFileTypes.includes("csv")}
              />
            </StyledInputCheckboxWrapper>
          </StyledSection>

          <StyledMessage>
            We regularly expand our single cell data corpus to improve results.
            Downloaded data and figures may differ in the future.
          </StyledMessage>

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
              data-test-id="dialog-download-button"
            >
              Download
            </DownloadButton>
          </StyledButtonContainer>
        </StyledModalContent>
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
function generateCsv({
  allChartProps,
  compare,
  selectedGenes,
  tissueName,
  availableFilters,
  selectedFilters,
  selectedOrganismId,
  availableOrganisms,
  selectedTissues,
}: {
  allChartProps: { [tissue: string]: ChartProps };
  compare: CompareId | undefined;
  selectedGenes: Props["selectedGenes"];
  tissueName: string;
  availableFilters: Partial<FilterDimensions>;
  selectedFilters: State["selectedFilters"];
  selectedOrganismId: string | null;
  availableOrganisms: OntologyTerm[] | null | undefined;
  selectedTissues: Props["selectedTissues"];
}) {
  const output: (string | number | undefined)[][] = [];

  // Create CSV comments and column names
  output.push(
    ...csvHeaders({
      compare,
      availableFilters,
      availableOrganisms,
      selectedFilters,
      selectedGenes,
      selectedOrganismId,
      selectedTissues,
    })
  );

  // Create a mapping of cell type IDs to a metadata array. (ex. "CL:00000" => [aggregated, female, male])
  const cellTypeIdToMetadataMapping = Object.values(
    buildCellTypeIdToMetadataMapping(tissueName, allChartProps)
  );

  for (const metadataForCellType of cellTypeIdToMetadataMapping) {
    for (const geneName of selectedGenes) {
      for (const metadata of metadataForCellType) {
        output.push(
          csvGeneExpressionRow({
            metadata,
            tissueName,
            allChartProps,
            geneName,
            compare,
          })
        );
      }
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
  compare,
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
  compare: CompareId | undefined;
  heatmapNode: HTMLDivElement;
  selectedTissues: Props["selectedTissues"];
  selectedGenes: Props["selectedGenes"];
  selectedCellTypes: Props["selectedCellTypes"];
  selectedFileTypes: ("png" | "svg" | "csv")[];
  availableFilters: Partial<FilterDimensions>;
  availableOrganisms: OntologyTerm[] | null | undefined;
  selectedOrganismId: string | null;
  selectedFilters: State["selectedFilters"];
  observer: MutationObserver;
  setDownloadStatus: Dispatch<
    SetStateAction<{
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
                  input = generateCsv({
                    allChartProps,
                    compare,
                    selectedGenes,
                    tissueName,
                    availableFilters,
                    selectedFilters,
                    selectedOrganismId,
                    availableOrganisms,
                    selectedTissues,
                  });
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
        file_type: selectedFileTypes,
        version:
          selectedFileTypes.length === 1 && selectedFileTypes[0] === "csv"
            ? "only_csv"
            : "includes_png_or_svg",

        dataset_filter: selectedFilters.datasets,
        disease_filter: selectedFilters.diseases,
        self_reported_ethnicity_filter: selectedFilters.ethnicities,
        sex_filter: selectedFilters.sexes,
        group_by_option: getCompareOptionNameById(compare),

        genes: selectedGenes,
        tissues: selectedTissues,
      });
    } catch (error) {
      console.error(error);
    }

    observer.disconnect();
    setEchartsRendererMode("canvas");
    setDownloadStatus({ isLoading: false });
  };
}
