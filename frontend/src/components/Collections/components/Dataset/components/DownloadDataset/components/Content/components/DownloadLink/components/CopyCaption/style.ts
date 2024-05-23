import styled from "@emotion/styled";
import { Link, fontBodyXs } from "@czi-sds/components";
import { primary400, textSecondary } from "src/common/theme";

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

export const StyledLink = styled(Link)`
  color: ${primary400};
`;
