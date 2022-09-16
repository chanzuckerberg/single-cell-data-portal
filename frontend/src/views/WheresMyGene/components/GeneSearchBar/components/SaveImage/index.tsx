import { Button as BPButton, Classes, Intent } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { FormControlLabel } from "@material-ui/core";
import { Icon, RadioButton } from "czifui";
import download from "downloadjs";
import { toPng, toSvg } from "html-to-image";
import { useCallback, useState } from "react";
import {
  Section,
  Title,
} from "src/components/Collections/components/Dataset/components/DownloadDataset/components/Content/components/common/style";
import Modal from "src/components/common/Modal";
import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import { Label } from "../../style";
import { ButtonWrapper, StyledIconButton } from "../QuickSelect/style";

export const EXCLUDE_IN_SCREENSHOT_CLASS_NAME = "screenshot-exclude";

const screenshotFilter = (domNode: HTMLElement): boolean => {
  return (
    !domNode.classList?.contains(EXCLUDE_IN_SCREENSHOT_CLASS_NAME) &&
    domNode.tagName !== "NOSCRIPT"
  );
};

const DownloadButton = styled(BPButton)`
  margin-right: ${PT_GRID_SIZE_PX}px;
`;

export default function SaveImage({
  selectedTissues,
  selectedGenes,
}: {
  selectedTissues: Array<string>;
  selectedGenes: Array<string>;
}): JSX.Element {
  const [isOpen, setIsOpen] = useState(false);
  const [fileType, setFileType] = useState<"png" | "svg">("png");
  const handleButtonClick = useCallback(() => {
    setIsOpen(!isOpen);
  }, [isOpen]);

  const saveImage = useCallback(async () => {
    try {
      const heatmapNode = document.getElementById("view") as HTMLCanvasElement;
      heatmapNode.classList.add("CLONED");
      const convertHTMLtoImage = fileType === "png" ? toPng : toSvg;
      const imageURL = await convertHTMLtoImage(heatmapNode, {
        backgroundColor: "white",
        clonedClassName: "CLONED",
        filter: screenshotFilter,
        style: {
          width: "min-content",
        },
      });
      heatmapNode.classList.remove("CLONED");
      download(
        imageURL,
        `CELLxGENE_gene_expression_${selectedTissues[0]}.${fileType}`
      );
    } catch (error) {
      console.error(error);
    }
  }, [fileType, selectedTissues]);

  console.log(selectedTissues, selectedGenes);

  return (
    <>
      <ButtonWrapper>
        <Label>Download</Label>
        <StyledIconButton
          disabled={selectedTissues.length !== 1 || selectedGenes.length === 0}
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
          <DownloadButton intent={Intent.PRIMARY} onClick={saveImage}>
            Download
          </DownloadButton>
        </div>
      </Modal>
    </>
  );
}
