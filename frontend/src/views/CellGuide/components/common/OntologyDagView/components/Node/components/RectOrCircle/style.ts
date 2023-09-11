import { CommonThemeProps } from "@czi-sds/components";
import styled from "@emotion/styled";
import { secondaryColor } from "../../../../common/constants";

interface NodeWrapperProps {
  columnGap: number;
}
export const NodeWrapper = styled.div<NodeWrapperProps>`
  display: flex;
  flex-direction: row;
  column-gap: ${(props) => `${props.columnGap}px`};
  position: relative;
`;

interface BorderProps extends CommonThemeProps {
  width: number;
  height: number;
}

export const Border = styled.div<BorderProps>`
  border-left: 1px solid ${secondaryColor};
  border-right: 1px solid ${secondaryColor};
  width: ${(props) => `${props.width}px`};
  height: ${(props) => `${props.height}px`};
`;
