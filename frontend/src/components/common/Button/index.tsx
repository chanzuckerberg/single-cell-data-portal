import styled from "@emotion/styled";
import {
  ButtonProps as RawButtonProps,
  CommonThemeProps,
  Button as RawButton,
  fontBodyS,
} from "@czi-sds/components";
import {
  grayWhite,
  spacesM,
  spacesXs,
  textPrimary,
  textSecondary,
} from "src/common/theme";

type ButtonProps = RawButtonProps & CommonThemeProps;

export const Button = styled(RawButton)`
  ${squarePrimary}
  ${minimalSecondary}
`;

function squarePrimary(props: ButtonProps) {
  const { sdsStyle, sdsType } = props;

  if (sdsStyle !== "square" || sdsType !== "primary") return;

  return `
    ${commonStyle()}
    min-width: 80px;
    color: ${grayWhite()};
    padding: ${spacesXs(props)}px ${spacesM(props)}px;
  `;
}

function minimalSecondary(props: ButtonProps) {
  const { sdsStyle, sdsType } = props;

  if (sdsStyle !== "minimal" || sdsType !== "secondary") return;

  return `
    ${commonStyle()}
    color: ${textSecondary(props)};

    &:hover {
      color: ${textPrimary(props)};
    }
  `;
}

function commonStyle() {
  return `
    ${fontBodyS}
    font-weight: 500;
  `;
}
