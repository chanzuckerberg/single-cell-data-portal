import { Button as BPButton, Classes, Intent } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { FormControlLabel } from "@material-ui/core";
import { Icon, RadioButton } from "czifui";
import { toPng, toSvg } from "html-to-image";
import { useCallback, useState } from "react";
import {
  Section,
  Title,
} from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/common/style";
import Modal from "src/components/common/Modal";
import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import { CellType } from "src/views/WheresMyGene/common/types";
import { getHeatmapHeight, getHeatmapWidth } from "../../../HeatMap/utils";
import { Label } from "../../style";
import { ButtonWrapper, StyledIconButton } from "../QuickSelect/style";

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

const DownloadButton = styled(BPButton)`
  margin-right: ${PT_GRID_SIZE_PX}px;
`;
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
}: {
  selectedTissues: Array<string>;
  selectedGenes: Array<string>;
  selectedCellTypes: { [tissue: string]: CellType[] };
}): JSX.Element {
  const [isOpen, setIsOpen] = useState(false);
  const [fileType, setFileType] = useState<"png" | "svg">("png");
  const handleButtonClick = useCallback(() => {
    setIsOpen(!isOpen);
  }, [isOpen]);

  const handleDownload = useCallback(async () => {
    try {
      const heatmapNode = document.getElementById("view") as HTMLCanvasElement;
      heatmapNode.classList.add("CLONED");
      const isPNG = fileType === "png";
      const convertHTMLtoImage = isPNG ? toPng : toSvg;
      const images = await Promise.all(
        selectedTissues.map(async (tissue) => {
          const imageURL = await convertHTMLtoImage(heatmapNode, {
            backgroundColor: "white",
            filter: screenshotFilter(tissue),
            height: getHeatmapHeight(selectedCellTypes[tissue]) + 200,
            pixelRatio: 1,
            width: getHeatmapWidth(selectedGenes) + 200,
          });

          return {
            input: isPNG
              ? base64URLToArrayBuffer(imageURL)
              : decodeURIComponent(imageURL.split(",")[1]),
            name: `${tissue}.${fileType}`,
          };
        })
      );
      heatmapNode.classList.remove("CLONED");
      const { downloadZip } = await import("client-zip");
      const blob = await downloadZip(images).blob();
      const link = document.createElement("a");
      link.href = URL.createObjectURL(blob);
      link.download = `CELLxGENE_1PX_gene_expression.zip`;
      link.click();
      link.remove();
    } catch (error) {
      console.error(error);
      // this could probably be handled better
      return { input: "", name: "" };
    }
  }, [fileType, selectedCellTypes, selectedTissues, selectedGenes]);

  return (
    <>
      <ButtonWrapper>
        <Label>Download</Label>
        <StyledIconButton
          disabled={selectedTissues.length === 0 || selectedGenes.length === 0}
          data-test-id={"download-button"}
          onClick={handleButtonClick}
          sdsType="primary"
          sdsSize="medium"
        >
          <Icon sdsIcon="download" sdsSize="l" sdsType="iconButton" />
        </StyledIconButton>
      </ButtonWrapper>
      <Modal
        isOpen={isOpen}
        title="Download Figure"
        onClose={handleButtonClick}
      >
        <div className={Classes.DIALOG_BODY}>
          <Section>
            <Title>Image Format</Title>
            <FormControlLabel
              control={
                <RadioButton
                  stage={fileType === "png" ? "checked" : "unchecked"}
                />
              }
              label=".png"
              onChange={() => setFileType("png")}
            />

            <FormControlLabel
              control={
                <RadioButton
                  stage={fileType === "svg" ? "checked" : "unchecked"}
                />
              }
              label=".svg"
              onChange={() => setFileType("svg")}
            />
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
