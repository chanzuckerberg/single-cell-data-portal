import React, { FC, useState } from "react";
import { COOKIE_KEYS } from "src/common/constants/cookieKeys";
import { getCookie } from "src/common/cookies/getCookies";
import { BOOLEAN, setCookie } from "src/common/cookies/setCookies";
import { ButtonWrapper, Link, NoButton, OKButton, Wrapper } from "./style";

const TOS_LINK = "https://cellxgene.cziscience.com/static/deploy/tos.html";
const PRIVACY_LINK =
  "https://cellxgene.cziscience.com/static/deploy/privacy.html";

const CookieBanner: FC = () => {
  const [isHidden, setIsHidden] = useState(getHasClickedOnBanner());

  if (isHidden) return null;

  function handleOKClick() {
    setCookie(COOKIE_KEYS.COOKIES_ACCEPTED, BOOLEAN.TRUE);
    setIsHidden(true);
  }

  function handleNoClick() {
    setCookie(COOKIE_KEYS.COOKIES_ACCEPTED, BOOLEAN.FALSE);
    setIsHidden(true);
  }

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
        <OKButton onClick={handleOKClick}>I&lsquo;m OK with cookies</OKButton>
        <NoButton onClick={handleNoClick}>No thanks</NoButton>
      </ButtonWrapper>
    </Wrapper>
  );
};

function getHasClickedOnBanner() {
  return getCookie(COOKIE_KEYS.COOKIES_ACCEPTED) !== "";
}

export default CookieBanner;
