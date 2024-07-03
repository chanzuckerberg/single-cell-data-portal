import styled from "@emotion/styled";
import { fontBodyXs, Tag } from "@czi-sds/components";
import { fontWeightSemibold, gray100, gray500 } from "src/common/theme";

export const CountAndTotal = styled(Tag)`
  &.MuiChip-root {
    background-color: ${gray100};
    margin: 0;

    .MuiChip-label {
      ${fontBodyXs}
      color: ${gray500};
      font-weight: ${fontWeightSemibold};
    }
  }
`;
