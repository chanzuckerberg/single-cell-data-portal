import React, { FC } from "react";
import { ButtonWrapper, Link, NoButton, OKButton, Wrapper } from "./style";

const TOS_LINK = "https://cellxgene.cziscience.com/static/deploy/tos.html";
const PRIVACY_LINK =
  "https://cellxgene.cziscience.com/static/deploy/privacy.html";

const CookieBanner: FC = () => {
  return (
    <Wrapper>
      By using this site, you are agreeing to our{" "}
      <Link href={TOS_LINK} target="_blank" rel="noreferrer">
        terms of service
      </Link>
      . We use cookies to help us improve our future efforts, and we also use
      necessary cookies to make our site work. To learn more, read our{" "}
      <Link href={PRIVACY_LINK} target="_blank" rel="noreferrer">
        privacy policy
      </Link>
      .
      <ButtonWrapper>
        <OKButton>I&lsquo;m OK with cookies</OKButton>
        <NoButton>No thanks</NoButton>
      </ButtonWrapper>
    </Wrapper>
  );
};

export default CookieBanner;
