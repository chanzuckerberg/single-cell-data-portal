import styled from "@emotion/styled";
import { fontBodyS, getColors } from "czifui";

export const Header = styled.h1`
  margin-bottom: 12px;
  font-size: 36px;
  font-weight: bold;
`;

export const Details = styled.p`
  ${fontBodyS}

  letter-spacing: 0;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;

export const Content = styled.div`
  display: flex;
  justify-content: space-between;
`;

export const Step3Details = styled.p`
  margin-bottom: 8px;
`;
