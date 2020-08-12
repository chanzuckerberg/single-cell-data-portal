import React, { FC, useState } from "react";
import { LOCAL_STORAGE_KEYS } from "src/common/constants/localStorageKeys";
import { get } from "src/common/localStorage/get";
import { BOOLEAN, set } from "src/common/localStorage/set";
import { ButtonWrapper, Link, NoButton, OKButton, Wrapper } from "./style";

const TOS_LINK = "https://cellxgene.cziscience.com/static/deploy/tos.html";
const PRIVACY_LINK =
  "https://cellxgene.cziscience.com/static/deploy/privacy.html";

const CookieBanner: FC = () => {
  const [isHidden, setIsHidden] = useState(hasClickedOnBanner());

  if (isHidden) return null;

  function handleOKClick() {
    set(LOCAL_STORAGE_KEYS.COOKIE_DECISION, BOOLEAN.TRUE);
    setIsHidden(true);
  }

  function handleNoClick() {
    set(LOCAL_STORAGE_KEYS.COOKIE_DECISION, BOOLEAN.FALSE);
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

function hasClickedOnBanner() {
  return get(LOCAL_STORAGE_KEYS.COOKIE_DECISION) !== null;
}

export default CookieBanner;
