import styled from "@emotion/styled";
import { CommonThemeProps, fontCapsXxs, fontHeaderL, getColors } from "czifui";

export const Wrapper = styled.div`
  /* Number is 40px wide, gap: 16px, Content: 200px */
  gap: 16px;
  height: 100%;

  border-radius: 4px;

  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

    return `
      background: ${colors?.gray[200]};

      border-radius: 5px;
    `;
  }}
`;

export const Header = styled.span`
  ${fontCapsXxs}

  display: block;

  padding: 24px 24px 4px 24px;

  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;

export const Details = styled.span`
  ${fontHeaderL}

  letter-spacing: 0;
  display: inline-block;
  padding-left: 24px;

  color: black;
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
