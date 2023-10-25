import styled from "@emotion/styled";
import { CommonThemeProps, getColors } from "@czi-sds/components";

interface DotProps extends CommonThemeProps {
  size: number;
}

export const Wrapper = styled.div`
  padding: 0 10px;
`;

export const Dot = styled.span`
  border-radius: 50%;

  ${(props: DotProps) => {
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
  margin-bottom: 3px;
  justify-content: space-between;
  align-items: center;
`;
