import styled from "@emotion/styled";
import { fontBodyS, getColors } from "@czi-sds/components";

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

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;
