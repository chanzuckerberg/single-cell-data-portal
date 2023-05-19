import styled from "@emotion/styled";
import { ButtonIcon, CommonThemeProps, getColors } from "czifui";

const primary600 = (props: CommonThemeProps) => getColors(props)?.primary[600];

export const RevisionButton = styled(ButtonIcon)`
  &:hover {
    color: ${primary600};
  }

  &:active {
    color: ${primary600};
  }

  .MuiSvgIcon-root {
    height: 16px;
    width: 16px;
  }
`;
