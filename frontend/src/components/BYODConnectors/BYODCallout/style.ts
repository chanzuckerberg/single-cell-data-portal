import styled from "@emotion/styled";
import { Callout, fontBodyS } from "@czi-sds/components";
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
