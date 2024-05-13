import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../SaveExport";
import { StyledButtonDiv, StyledButtonIcon, StyledLabel } from "./style";
import { useConnect } from "./connect";

export default function ShareButton(): JSX.Element {
  const { copyShareUrl, selectedGenes, isDisabled } = useConnect();

  return (
    <StyledButtonDiv className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
      <StyledLabel>Share</StyledLabel>
      <StyledButtonIcon
        data-testid={"share-button"}
        onClick={copyShareUrl}
        sdsSize="medium"
        sdsType="primary"
        sdsIcon="share"
        disabled={selectedGenes.length === 0 || isDisabled}
      />
    </StyledButtonDiv>
  );
}
