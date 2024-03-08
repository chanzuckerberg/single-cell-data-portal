import styled from "@emotion/styled";
import {
  ButtonIcon as SDSButtonIcon,
  CommonThemeProps,
} from "@czi-sds/components";
import { cornersL, spacesL, spacesXxs, textPrimary } from "src/common/theme";
import { css } from "@emotion/react";

interface ButtonIconProps extends CommonThemeProps {
  open: boolean;
}

export const ButtonIcon = styled(SDSButtonIcon)<ButtonIconProps>`
  border-radius: ${cornersL}px;
  padding: ${spacesXxs}px;

  &:hover {
    color: ${textPrimary};
  }

  .MuiSvgIcon-root {
    height: ${spacesL}px;
    width: ${spacesL}px;
  }

  ${(props) =>
    props.open &&
    css`
      color: ${textPrimary(props)};
    `}
`;
