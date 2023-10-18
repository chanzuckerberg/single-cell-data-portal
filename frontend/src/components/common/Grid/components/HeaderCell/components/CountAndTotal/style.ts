import styled from "@emotion/styled";
import { fontBodyXs, Tag } from "@czi-sds/components";
import { gray100, gray600 } from "src/common/theme";

export const CountAndTotal = styled(Tag)`
  &.MuiChip-root {
    background-color: ${gray100};
    margin: 0;

    .MuiChip-label {
      ${fontBodyXs}
      color: ${gray600};
      font-weight: 500;
      letter-spacing: -0.003em;
    }
  }
`;
