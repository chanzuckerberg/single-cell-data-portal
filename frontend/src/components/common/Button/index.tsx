import styled from "@emotion/styled";
import {
  Button as RawButton,
  ButtonProps as RawButtonProps,
  CommonThemeProps,
  fontBodyS,
} from "@czi-sds/components";
import {
  grayWhite,
  spacesM,
  spacesXs,
  textPrimary,
  textSecondary,
} from "src/common/theme";
import { css, SerializedStyles } from "@emotion/react";

type ButtonProps = RawButtonProps & CommonThemeProps;

export const Button = styled(RawButton)`
  ${squarePrimary}
  ${minimalSecondary}
`;

function squarePrimary(props: ButtonProps): SerializedStyles | undefined {
  const { sdsStyle, sdsType } = props;

  if (sdsStyle !== "square" || sdsType !== "primary") return;

  return css`
    ${commonStyle(props)}
    min-width: 80px;
    color: ${grayWhite()};
    padding: ${spacesXs(props)}px ${spacesM(props)}px;
  `;
}

function minimalSecondary(props: ButtonProps): SerializedStyles | undefined {
  const { sdsStyle, sdsType } = props;

  if (sdsStyle !== "minimal" || sdsType !== "secondary") return;

  return css`
    ${commonStyle(props)}
    color: ${textSecondary(props)};

    &:hover {
      color: ${textPrimary(props)};
    }
  `;
}

function commonStyle(props: ButtonProps): SerializedStyles {
  return css`
    ${fontBodyS(props)}
    font-weight: 500;
  `;
}
