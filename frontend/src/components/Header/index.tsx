import React, { FC } from "react";
import { HomepageLink } from "../common/HomepageLink";
import { Wrapper } from "./style";

const Header: FC = () => {
  return (
    <Wrapper>
      <HomepageLink />
    </Wrapper>
  );
};

export default Header;
