import styled from "@emotion/styled";
import { Callout, Button, fontBodyS } from "@czi-sds/components";
import {
  fontWeightSemibold,
  fontWeightRegular,
  spacesXxs,
  spacesM,
  primary100,
} from "src/common/theme";

export const StyledCallout = styled(Callout)`
  width: 100%;
  margin-right: ${spacesM}px;
  box-sizing: border-box;

  /* Ensure the callout doesn't expand beyond its container */
  max-width: calc(100% - ${spacesM}px);

  /* Fix background color for info intent */
  background-color: ${primary100};

  /* Target the text content */
  > div {
    font-weight: ${fontWeightRegular};
  }
`;

export const CalloutTitle = styled.span`
  ${fontBodyS}
  font-weight: ${fontWeightSemibold};
`;

export const CalloutTextWrapper = styled.div`
  margin-bottom: ${spacesXxs}px;
`;

// We need this because this version of SDS does not respect the isAllCaps=false prop.
export const StyledButton = styled(Button)`
  text-transform: none;

  /* Apply SDS body-s-500 typography */
  font-size: 14px;
  font-weight: 500;
  letter-spacing: 0px;
  line-height: 24px;
  font-family:
    var(--font-inter),
    Inter,
    -apple-system,
    BlinkMacSystemFont,
    Segoe UI,
    Roboto,
    Helvetica Neue,
    Helvetica,
    Arial,
    sans-serif;
`;
