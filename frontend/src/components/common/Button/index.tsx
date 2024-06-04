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
  success400,
  success500,
  success600,
  textPrimary,
  textSecondary,
} from "src/common/theme";
import { css, SerializedStyles } from "@emotion/react";

type ButtonProps = RawButtonProps & CommonThemeProps;

export const Button = styled(RawButton)<ButtonProps>`
  ${squarePrimary}
  ${squareSecondary}
  ${minimalSecondary}
  ${minimalPrimary}
`;

function squarePrimary(props: ButtonProps): SerializedStyles | undefined {
  const { color, sdsStyle, sdsType } = props;

  if (sdsStyle !== "square" || sdsType !== "primary") return;

  if (color === "success") {
    return css`
      ${commonStyle(props)}
      ${squarePrimarySuccessStyle(props)}
    `;
  }

  return css`
    ${commonStyle(props)}
    min-width: 80px;
    color: ${grayWhite()};
    padding: ${spacesXs(props)}px ${spacesM(props)}px;
  `;
}

function squareSecondary(props: ButtonProps): SerializedStyles | undefined {
  const { sdsStyle, sdsType } = props;

  if (sdsStyle !== "square" || sdsType !== "secondary") return;

  return css`
    ${commonStyle(props)}
    min-width: 80px;
    padding: ${spacesXs(props)}px ${spacesM(props)}px;
  `;
}

function minimalPrimary(props: ButtonProps): SerializedStyles | undefined {
  const { sdsStyle, sdsType } = props;

  if (sdsStyle !== "minimal" || sdsType !== "primary") return;

  return css`
    ${commonStyle(props)}
    ${textTransformStyle(props)}
  `;
}

function minimalSecondary(props: ButtonProps): SerializedStyles | undefined {
  const { sdsStyle, sdsType } = props;

  if (sdsStyle !== "minimal" || sdsType !== "secondary") return;

  return css`
    ${commonStyle(props)}
    ${textTransformStyle(props)}
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

function squarePrimarySuccessStyle(props: ButtonProps): SerializedStyles {
  return css`
    box-shadow: inset 0 0 0 1px ${success400(props)};

    &:hover {
      background-color: ${success500(props)};
      box-shadow: inset 0 0 0 1px ${success500(props)};
    }

    &:active {
      background-color: ${success600(props)};
      box-shadow: inset 0 0 0 1px ${success600(props)};
    }
  `;
}

// Custom fix for SDS button sdsType "minimal" button; isAllCaps is expected to be an optional style see https://github.com/chanzuckerberg/sci-components/blob/main/packages/components/src/core/Button/index.tsx#L68 and https://github.com/chanzuckerberg/sci-components/blob/main/packages/components/src/core/Button/style.ts#L219.
function textTransformStyle(props: ButtonProps): SerializedStyles | undefined {
  const { isAllCaps } = props;
  if (isAllCaps || isAllCaps === undefined) return;
  return css`
    text-transform: none;
  `;
}
