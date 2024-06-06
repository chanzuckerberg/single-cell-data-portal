import { StyledButtonDiv, StyledLabel, StyledButtonIcon } from "./style";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../SaveExport";
import { useConnect } from "./connect";
import { CITATION_BUTTON_LABEL } from "./constants";

export default function CitationButton(): JSX.Element {
  const { copyCitation, isDisabled } = useConnect();

  return (
    <StyledButtonDiv className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
      <StyledLabel>{CITATION_BUTTON_LABEL}</StyledLabel>
      <StyledButtonIcon
        data-testid={"citation-button"}
        onClick={copyCitation}
        sdsSize="medium"
        sdsType="primary"
        icon="Quote"
        disabled={isDisabled}
      />
    </StyledButtonDiv>
  );
}
