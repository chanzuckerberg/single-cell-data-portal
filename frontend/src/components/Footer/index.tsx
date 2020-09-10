import React, { FC } from "react";
import { MainWrapper, StyledAnchor, Wrapper } from "./style";

const Footer: FC = () => {
  return (
    <Wrapper>
      <MainWrapper>
        © {new Date().getFullYear()} cellxgene
        <StyledAnchor
          href={"mailto:cellxgene@chanzuckerberg.com"}
          target="_blank"
          rel="noopener"
        >
          Contact
        </StyledAnchor>
        <StyledAnchor
          href={
            "https://cellxgene.cziscience.com/static/cellxgene/deploy/tos.html"
          }
          target="_blank"
          rel="noopener"
        >
          Terms of Service
        </StyledAnchor>
        <StyledAnchor
          href={
            "https://cellxgene.cziscience.com/static/cellxgene/deploy/privacy.html"
          }
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
