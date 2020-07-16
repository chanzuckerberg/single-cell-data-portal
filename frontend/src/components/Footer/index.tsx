import React, { FC } from "react";
import { HomepageLink } from "../common/HomepageLink";
import { Wrapper } from "./style";

const Footer: FC = () => {
  return (
    <Wrapper>
      <HomepageLink small />
    </Wrapper>
  );
};

export default Footer;
