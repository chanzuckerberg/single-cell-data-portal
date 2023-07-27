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
  X_AXIS_CHART_HEIGHT_PX,
  X_AXIS_HOVER_CONTAINER_HEIGHT_PX,
} from "src/views/WheresMyGene/common/constants";
import {
  CellType,
  ChartProps,
  Tissue,
} from "src/views/WheresMyGene/common/types";

import {
  getHeatmapHeight,
  getHeatmapWidth,
  hyphenize,
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
  DOWNLOAD_MODAL_WIDTH_PX,
  DOWNLOAD_MODAL_PADDING,
} from "./style";
import {
  applyAttributes,
  NAME_SPACE_URI,
  renderDots,
  renderLegend,
  renderXAxis,
  renderYAxis,
} from "./svgUtils";
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
  getCurrentDate,
} from "./csvUtils";
import { InputCheckbox } from "@czi-sds/components";
import {
  DATA_MESSAGE_BANNER_HEIGHT_PX,
  DATA_MESSAGE_BANNER_ID,
  DATA_MESSAGE_BANNER_WIDTH_PX,
  UnderlyingDataChangeBanner,
} from "./ExportBanner";
import {
  CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX,
  CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX,
} from "src/components/Layout/style";
import { LEGEND_HEIGHT_PX } from "../../../InfoPanel/components/Legend/style";
import { LEGEND_MARGIN_BOTTOM_PX } from "src/views/WheresMyGene/style";
import { CHART_PADDING_PX } from "../../../HeatMap/style";

let heatmapContainerScrollTop: number | undefined;

const MUTATION_OBSERVER_TIMEOUT = 3 * 1000;
const CLONED_CLASS = "CLONED";

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

    // We only need the observer to listen to changes to the heatmap
    // If the observer listens on the whole view then it will detect the banner injection and trigger download twice
    const heatmapNode = document.getElementById(
      HEATMAP_CONTAINER_ID
    ) as HTMLDivElement;

    // This is the node that will be used to generate the exports
    const viewNode = document.getElementById("view") as HTMLDivElement;

    // Options for the observer (which mutations to observe)
    const config: MutationObserverInit = {
      childList: true,
      subtree: true,
      attributes: true, // triggers the the observer if adding a class, useful for forcing the observer to trigger
    };

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
        heatmapNode: viewNode,
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

    // This will make a change to the heatmap dom which triggers the observer to start the download
    if (
      selectedFileTypes.includes("svg") ||
      selectedFileTypes.includes("png")
    ) {
      setEchartsRendererMode("svg");
    } else {
      // Kind of a hack to modify the DOM to trigger the observer
      heatmapNode.classList.add("is-downloading");
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
          data-testid="download-button"
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
                caption="Image files that can be used on websites, have transparent backgrounds, and are suitable for graphics and digital images."
                onChange={() => selectFileType("png")}
                checked={selectedFileTypes.includes("png")}
                data-testid="png-checkbox"
              />
            </StyledInputCheckboxWrapper>
            <StyledInputCheckboxWrapper>
              <InputCheckbox
                label="SVG"
                caption="Vector image files that can be scaled without losing quality and are most suitable for scientific publications, presentations, and data visualization."
                onChange={() => selectFileType("svg")}
                checked={selectedFileTypes.includes("svg")}
                data-testid="svg-checkbox"
              />
            </StyledInputCheckboxWrapper>
          </StyledSection>

          <StyledSection>
            <StyledTitle>DATA</StyledTitle>
            <StyledInputCheckboxWrapper>
              <InputCheckbox
                label="CSV"
                caption="Plain text files that store tabular data by separating values with commas, making them easy to read and manipulate with spreadsheet software."
                onChange={() => selectFileType("csv")}
                checked={selectedFileTypes.includes("csv")}
                data-testid="csv-checkbox"
              />
            </StyledInputCheckboxWrapper>
          </StyledSection>

          <StyledMessage>
            <UnderlyingDataChangeBanner
              centered={false} // Decided by UX team that text should only be centered for exports
              width={DOWNLOAD_MODAL_WIDTH_PX - DOWNLOAD_MODAL_PADDING * 2} // setting the width of the banner svg to the modal's content width
            />
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
              data-testid="dialog-download-button"
            >
              Download
            </DownloadButton>
          </StyledButtonContainer>
        </StyledModalContent>
      </StyledModal>
    </>
  );
}

function generateSvg({
  svg,
  tissueName,
  heatmapHeight,
  heatmapWidth,
}: {
  svg: string;
  tissueName: Tissue;
  heatmapHeight: number;
  heatmapWidth: number;
}) {
  const heatmapNode = new DOMParser().parseFromString(svg, "image/svg+xml");
  const heatmapContainer = heatmapNode
    .querySelector("foreignObject")
    ?.querySelector("div");

  // This is an svg element created in dom when downloadStatus.isLoading is true
  const banner = document.getElementById(DATA_MESSAGE_BANNER_ID);

  if (!banner) return "";

  // Used to create room for the banner
  const paddedBannerHeight =
    DATA_MESSAGE_BANNER_HEIGHT_PX + CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX;
  const paddedBannerWidth =
    DATA_MESSAGE_BANNER_WIDTH_PX + CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX * 2;

  // Render elements to SVG
  const xAxisSvg = renderXAxis({
    heatmapContainer,
    yOffset: paddedBannerHeight,
  });
  const yAxisSvg = renderYAxis({
    heatmapContainer,
    heatmapHeight,
    tissueName,
    yOffset: paddedBannerHeight,
  });
  const dotsSvg = renderDots({
    heatmapContainer,
    tissueName,
    yOffset: paddedBannerHeight,
  });
  const legendSvg = renderLegend({
    heatmapContainer,
    yOffset: paddedBannerHeight,
  });

  const svgWidth =
    heatmapWidth +
    Y_AXIS_CHART_WIDTH_PX +
    CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX * 2;

  const finalSvg = document.createElementNS(NAME_SPACE_URI, "svg");
  applyAttributes(finalSvg, {
    width: svgWidth < paddedBannerWidth ? paddedBannerWidth : svgWidth, // Use the banner width as the minimum final svg width
    height:
      heatmapHeight +
      paddedBannerHeight +
      X_AXIS_CHART_HEIGHT_PX +
      CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX * 2,
  });

  // Center banner using full svg width
  applyAttributes(banner, {
    x:
      svgWidth > paddedBannerWidth
        ? svgWidth / 2 - DATA_MESSAGE_BANNER_WIDTH_PX / 2
        : CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX,
    y: CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX,
  });

  // Required for valid SVG
  finalSvg.setAttributeNS(
    "http://www.w3.org/2000/xmlns/",
    "xmlns",
    NAME_SPACE_URI
  );

  // Append elements to final SVG
  finalSvg.append(banner.cloneNode(true) || ""); // We need to clone this or else the element itself will be moved and not be accessible for subsequent tissues
  finalSvg.append(legendSvg || "");
  finalSvg.append(xAxisSvg || "");
  finalSvg.append(yAxisSvg || "");
  finalSvg.append(dotsSvg || "");

  return finalSvg.outerHTML;
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
  availableFilters,
  selectedFilters,
  selectedOrganismId,
  availableOrganisms,
  selectedTissues,
}: {
  allChartProps: { [tissue: string]: ChartProps };
  compare: CompareId | undefined;
  selectedGenes: Props["selectedGenes"];
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

  for (const tissueName of selectedTissues) {
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
  }

  return csvStringify(output);
}

async function generateImage({
  fileType,
  heatmapNode,
  heatmapHeight,
  heatmapWidth,
  tissueName,
  isMultipleTissues = false,
}: {
  fileType: string;
  heatmapNode: HTMLDivElement;
  heatmapHeight: number;
  heatmapWidth: number;
  tissueName: string;
  isMultipleTissues: boolean;
}): Promise<string | ArrayBuffer> {
  const convertHTMLtoImage = fileType === "png" ? toPng : toSvg;

  const imageURL = await convertHTMLtoImage(heatmapNode, {
    backgroundColor: "white",
    filter: screenshotFilter(tissueName),
    height:
      heatmapHeight +
      X_AXIS_HOVER_CONTAINER_HEIGHT_PX +
      X_AXIS_CHART_HEIGHT_PX +
      DATA_MESSAGE_BANNER_HEIGHT_PX +
      LEGEND_HEIGHT_PX +
      LEGEND_MARGIN_BOTTOM_PX +
      CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX * 2,
    pixelRatio: 4,
    width: heatmapNode.offsetWidth,
  });

  let input: string | ArrayBuffer = imageURL;

  if (fileType === "svg") {
    input = generateSvg({
      heatmapHeight,
      heatmapWidth,
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

      const heatmapWidth = getHeatmapWidth(selectedGenes);

      if (isPng) {
        // Add classes that are required for styling PNG
        // Adding this class to the heatmap causes the y-axis scrolling to jump but is required for image download
        heatmapNode.classList.add(CLONED_CLASS);
        document.getElementById("top-legend")?.classList.add(CLONED_CLASS);

        heatmapNode.style.width = `${
          heatmapWidth +
          Y_AXIS_CHART_WIDTH_PX +
          CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX * 2 +
          CHART_PADDING_PX * 2
        }px`;
      }

      const exports: ExportData[] = [];

      await Promise.all(
        selectedTissues.map(async (tissueName) => {
          // Handles if whitespace is in the tissue name for the element ID
          const formattedTissueName = hyphenize(tissueName);

          // Generate exports for images only
          return await Promise.all(
            selectedFileTypes.map(async (fileType) => {
              if (fileType === "png" || fileType === "svg") {
                const heatmapHeight = getHeatmapHeight(
                  selectedCellTypes[tissueName]
                );

                exports.push({
                  input: await generateImage({
                    fileType,
                    heatmapNode,
                    heatmapHeight,
                    heatmapWidth,
                    tissueName: formattedTissueName,
                    isMultipleTissues:
                      selectedTissues.length > 1 ||
                      selectedFileTypes.length > 1,
                  }),
                  name: `${formattedTissueName}.${fileType}`,
                });
              }
            })
          );
        })
      );

      // CSV has all tissues in one file
      if (selectedFileTypes.includes("csv")) {
        exports.push({
          input: generateCsv({
            allChartProps,
            compare,
            selectedGenes,
            availableFilters,
            selectedFilters,
            selectedOrganismId,
            availableOrganisms,
            selectedTissues,
          }),
          name:
            selectedTissues.length === 1 // If only one tissue is selected, use tissue name as filename
              ? `${selectedTissues[0]}.csv`
              : `CELLxGENE_gene_expression_${getCurrentDate()}.csv`,
        });
      }

      if (isPng) {
        // Remove classes that were required for styling PNG
        heatmapNode.classList.remove(CLONED_CLASS);
        document.getElementById("top-legend")?.classList.remove(CLONED_CLASS);

        //(thuang): #3569 Restore scrollTop position
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
        publication_filter: selectedFilters.publications,
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
