import React, { FC } from "react";
import logo from "src/common/images/logo.svg";
import { Logo as StyledLogo, LogoWrapper } from "./style";

interface Props {
  small?: boolean;
}

export const Logo: FC<Props> = ({ small }) => {
  return (
    <LogoWrapper>
      <StyledLogo size={small ? 100 : undefined} src={String(logo)} />
    </LogoWrapper>
  );
};
