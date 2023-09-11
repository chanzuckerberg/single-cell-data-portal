import { CommonThemeProps } from "@czi-sds/components";
import styled from "@emotion/styled";
import { secondaryColor } from "../../../../common/constants";

export const NodeWrapper = styled.div`
  display: flex;
  flex-direction: row;
  position: relative;
`;

interface BorderProps extends CommonThemeProps {
  width: number;
  height: number;
}

export const Border = styled.div<BorderProps>`
  border-left: 1px solid ${secondaryColor};
  background-color: #f8f8f8;
  width: ${(props) => `${props.width}px`};
  height: ${(props) => `${props.height / 2}px`};
  margin: auto 0;
`;
