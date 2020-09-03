import React, { FC } from "react";
import { HomepageLink } from "../common/HomepageLink";
import { MainWrapper, Wrapper } from "./style";

const Header: FC = () => {
  return (
    <Wrapper>
      <MainWrapper>
        <HomepageLink />
      </MainWrapper>
    </Wrapper>
  );
};

export default Header;
