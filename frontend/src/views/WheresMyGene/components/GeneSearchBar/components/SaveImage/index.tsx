import { Classes, Intent } from "@blueprintjs/core";
import { FormControlLabel } from "@mui/material";
import { InputRadio } from "czifui";
import { toPng, toSvg } from "html-to-image";
import { useCallback, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  Section,
  Title,
} from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/common/style";
import Modal from "src/components/common/Modal";
import { HEATMAP_CONTAINER_ID } from "src/views/WheresMyGene/common/constants";
import { CellType } from "src/views/WheresMyGene/common/types";
import {
  getHeatmapHeight,
  getHeatmapWidth,
  X_AXIS_CHART_HEIGHT_PX,
  Y_AXIS_CHART_WIDTH_PX,
} from "../../../HeatMap/utils";
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

export default function SaveImage({
  selectedTissues,
  selectedGenes,
  selectedCellTypes,
  setIsDownloading,
}: {
  selectedTissues: Array<string>;
  selectedGenes: Array<string>;
  selectedCellTypes: { [tissue: string]: CellType[] };
  setIsDownloading: (isDownloading: boolean) => void;
}): JSX.Element {
  const [isOpen, setIsOpen] = useState(false);
  const [fileType, setFileType] = useState<"png" | "svg">("png");
  const handleButtonClick = useCallback(() => {
    if (!isOpen) track(EVENTS.WMG_DOWNLOAD_CLICKED);
    setIsOpen(!isOpen);
  }, [isOpen]);

  const handleDownload = useCallback(async () => {
    setIsDownloading(true);
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
        selectedTissues.map(async (tissue) => {
          const imageURL = await convertHTMLtoImage(heatmapNode, {
            backgroundColor: "white",
            filter: screenshotFilter(tissue),
            height:
              getHeatmapHeight(selectedCellTypes[tissue]) +
              X_AXIS_CHART_HEIGHT_PX +
              120,
            pixelRatio: 4,
            width: heatmapNode.width,
          });
          // raw URI if only one tissue is selected
          const input =
            selectedTissues.length === 1
              ? imageURL
              : isPNG // otherwise, convert to array buffer if PNG
              ? base64URLToArrayBuffer(imageURL)
              : decodeURIComponent(imageURL.split(",")[1]);

          return {
            input,
            name: `${tissue}.${fileType}`,
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
        link.href = images[0].input as string;
        link.download = images[0].name;
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
  }, [fileType, selectedCellTypes, selectedTissues, selectedGenes]);

  return (
    <>
      <ButtonWrapper className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
        <Label>Download</Label>
        <StyledButtonIcon
          data-test-id="download-button"
          // TODO: put handleButtonClick when svgs are fixed
          onClick={handleDownload}
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
