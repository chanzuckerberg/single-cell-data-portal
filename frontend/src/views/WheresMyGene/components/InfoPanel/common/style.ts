import styled from "@emotion/styled";
import { fontBodyXxs, fontHeaderXxs, getColors } from "@czi-sds/components";

export const LowHigh = styled.div`
  display: flex;
  justify-content: space-between;

  /* Fixes bug where numbers were wrapping in PNG export */
  white-space: nowrap;
`;

export const Header = styled.h5`
  ${fontHeaderXxs}

  margin-bottom: 10px;
`;

export const Label = styled.label`
  ${fontBodyXxs}

  white-space: nowrap;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]}
    `;
  }}
`;

export const Content = styled.div`
  ${fontBodyXxs}

  width: 120px;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;
