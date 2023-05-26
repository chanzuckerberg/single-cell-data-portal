import styled from "@emotion/styled";
import { CommonThemeProps, getSpaces } from "@czi-sds/components";

import { Wrapper as RawWrapper } from "../common/style";

export const ViewOptionsWrapper = styled(RawWrapper)`
  ${(props: CommonThemeProps) => {
    const spaces = getSpaces(props);

    return `
      gap: ${spaces?.m}px;
    `;
  }}
`;
