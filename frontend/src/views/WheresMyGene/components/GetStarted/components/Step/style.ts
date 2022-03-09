import styled from "@emotion/styled";
import { fontBodyS, fontHeaderS, getColors } from "czifui";

export const Wrapper = styled.div`
  display: flex;
  /* Number is 40px wide, gap: 16px, Content: 200px */
  width: 256px;
  gap: 16px;
`;

export const Content = styled.div`
  width: 200px;
`;

export const Header = styled.h6`
  ${fontHeaderS}

  margin-bottom: 8px;
`;

export const Details = styled.p`
  ${fontBodyS}

  letter-spacing: 0;
  line-height: 18px;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;

export const Number = styled.span`
  border-radius: 50%;
  width: 40px;
  height: 40px;
  background-color: black;
  display: flex;
  justify-content: center;
  align-items: center;
`;

export const NumberContent = styled.span`
  font-size: 24px;
  font-weight: bold;
  color: white;
`;
