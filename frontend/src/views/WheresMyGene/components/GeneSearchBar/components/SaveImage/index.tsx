import { Classes, Intent } from "@blueprintjs/core";
import { FormControlLabel } from "@mui/material";
import { InputRadio } from "czifui";
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
import { CellType, Tissue } from "src/views/WheresMyGene/common/types";

import {
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
        fileType,
        heatmapNode,
        observer,
        selectedCellTypes,
        selectedGenes,
        selectedTissues,
        setEchartsRendererMode,
        setIsDownloading,
      }),
      MUTATION_OBSERVER_TIMEOUT
    );

    observer.observe(heatmapNode, config);

    setEchartsRendererMode("svg");
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

function download_({
  heatmapNode,
  selectedTissues,
  selectedGenes,
  selectedCellTypes,
  fileType,
  observer,
  setIsDownloading,
  setEchartsRendererMode,
}: {
  heatmapNode: HTMLDivElement;
  selectedTissues: Props["selectedTissues"];
  selectedGenes: Props["selectedGenes"];
  selectedCellTypes: Props["selectedCellTypes"];
  fileType: "png" | "svg";
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
            width: heatmapNode.offsetWidth,
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

    observer.disconnect();
    setIsDownloading(false);
    setEchartsRendererMode("canvas");
  };
}
