import styled from "@emotion/styled";
/* eslint-disable @typescript-eslint/ban-ts-comment --- TODO revisit lint error below */
// @ts-ignore - TODO revisit lint error
import { CommonThemeProps, getColors } from "czifui";

interface DotProps extends CommonThemeProps {
  size: number;
}

export const Wrapper = styled.div`
  margin-bottom: 24px;
`;

export const Dot = styled.span`
  border-radius: 50%;

  ${(props: DotProps) => {
    /* eslint-disable @typescript-eslint/ban-ts-comment --- TODO revisit lint error below */
    // @ts-ignore - TODO revisit lint error
    const colors = getColors(props);
    const { size } = props;

    return `
      background-color: ${colors?.gray[300]};
      width: ${size}px;
      height: ${size}px;
    `;
  }}
`;

export const Dots = styled.div`
  display: flex;
  margin-bottom: 10px;
  justify-content: space-between;
  align-items: center;
`;
