import React, { FC } from "react";
import { MainWrapper, StyledAnchor, Wrapper } from "./style";

const Footer: FC = () => {
  return (
    <Wrapper>
      <MainWrapper>
        © 2020 cellxgene
        <StyledAnchor
          key={"mailto:cellxgene@chanzuckerberg.com"}
          href={"mailto:cellxgene@chanzuckerberg.com"}
          target="_blank"
          rel="noopener"
        >
          Contact
        </StyledAnchor>
        <StyledAnchor
          key={"https://cellxgene.cziscience.com/static/cellxgene/deploy/tos.html"}
          href={"https://cellxgene.cziscience.com/static/cellxgene/deploy/tos.html"}
          target="_blank"
          rel="noopener"
        >
          Terms of Service
        </StyledAnchor>
        <StyledAnchor
          key={"https://cellxgene.cziscience.com/static/cellxgene/deploy/privacy.html"}
          href={"https://cellxgene.cziscience.com/static/cellxgene/deploy/privacy.html"}
          target="_blank"
          rel="noopener"
        >
          Privacy Policy
        </StyledAnchor>
      </MainWrapper>
    </Wrapper>
  );
};

export default Footer;
