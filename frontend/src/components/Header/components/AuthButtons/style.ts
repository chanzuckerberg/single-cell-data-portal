import styled from "@emotion/styled";
import { fontBodyS } from "@czi-sds/components";
import { gray500 } from "src/common/theme";

export const LogOutAnchor = styled.a`
  &:hover {
    text-decoration: none;
  }
`;

export const LogOutText = styled.div`
  ${fontBodyS}
`;

export const LogOutEmail = styled.div`
  ${fontBodyS}

  color: ${gray500};
`;
