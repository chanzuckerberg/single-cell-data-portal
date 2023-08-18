import { Classes, InputGroup } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { Callout } from "@czi-sds/components";
import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";

export const StyledInputGroup = styled(InputGroup)`
  & .${Classes.INPUT_ACTION} {
    top: 0;
  }
`;
export const FullWidthCallout = styled(Callout)`
  width: 100% !important;
  margin: ${PT_GRID_SIZE_PX * 2}px 0;
`;

export const APIDisclaimerP = styled.p`
  margin-top: 8px;
  width: 560px;
  color: ${GRAY.A};
`;
