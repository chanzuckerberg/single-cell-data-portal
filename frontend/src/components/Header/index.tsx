import React, { FC } from "react";
import { HomepageLink } from "../common/HomepageLink";
import AuthButtons from "./components/AuthButtons";
import { MainWrapper, Wrapper } from "./style";

const Header: FC = () => {
  return (
    <Wrapper>
      <MainWrapper>
        <HomepageLink />
        <AuthButtons />
      </MainWrapper>
    </Wrapper>
  );
};

export default Header;
