import { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
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
        <StyledAnchor href={ROUTES.TOS} target="_blank" rel="noopener">
          Terms of Service
        </StyledAnchor>
        <StyledAnchor href={ROUTES.PRIVACY} target="_blank" rel="noopener">
          Privacy Policy
        </StyledAnchor>
      </MainWrapper>
    </Wrapper>
  );
};

export default Footer;
