import styled from "@emotion/styled";
import { Button as SDSButton, CommonThemeProps } from "@czi-sds/components";
import { textPrimary } from "src/common/theme";
import { css } from "@emotion/react";

interface ButtonProps extends CommonThemeProps {
  open: boolean;
}

export const Button = styled(SDSButton, {
  shouldForwardProp: (prop) => prop !== "open",
})<ButtonProps>`
  ${(props) =>
    props.open &&
    css`
      .MuiSvgIcon-root {
        color: ${textPrimary(props)};
      }
    `}
`;
