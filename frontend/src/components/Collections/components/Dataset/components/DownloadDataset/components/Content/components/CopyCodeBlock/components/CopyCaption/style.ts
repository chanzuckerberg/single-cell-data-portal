import styled from "@emotion/styled";
import { fontBodyXs } from "@czi-sds/components";
import { textSecondary } from "src/common/theme";

export const Caption = styled.div`
  ${fontBodyXs}
  color: ${textSecondary};

  p {
    margin: 0;

    & + p {
      margin-top: 20px;
    }
  }

  b {
    font-weight: 500;
  }
`;
