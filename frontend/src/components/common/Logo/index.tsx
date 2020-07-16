import React, { FC } from "react";
import { CorporaWrapper, Logo as StyledLogo, LogoWrapper } from "./style";

interface Props {
  small?: boolean;
}

export const Logo: FC<Props> = ({ small }) => {
  return (
    <LogoWrapper>
      <StyledLogo
        size={small ? 20 : undefined}
        src="https://chanzuckerberg.com/wp-content/themes/czi/img/logo-minified.svg"
      />
      <CorporaWrapper small={small}>Corpora</CorporaWrapper>
    </LogoWrapper>
  );
};
