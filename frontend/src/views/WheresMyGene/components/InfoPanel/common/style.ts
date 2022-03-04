import styled from "@emotion/styled";
import { fontBodyXxs, fontHeaderXxs, getColors } from "czifui";

export const LowHigh = styled.div`
  display: flex;
  justify-content: space-between;
`;

export const Header = styled.h5`
  ${fontHeaderXxs}

  margin-bottom: 10px;
`;

export const Content = styled.div`
  ${fontBodyXxs}

  width: 100px;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;
