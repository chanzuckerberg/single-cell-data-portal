import Image from "next/image";
import Script from "next/script";
import { memo } from "react";
import {
  BOTTOM_BANNER_ID,
  NewsletterModal,
  StyledBanner,
  StyledBottomBannerWrapper,
  StyledLink,
  HeaderContainer,
  HiddenHubSpotForm,
  FooterContentWrapper,
  StyledCloseButtonIcon,
} from "./style";
import cellxgeneLogoSvg from "src/common/images/CellxGene.svg";
import Head from "next/head";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SaveExport";
import { noop } from "src/common/constants/utils";
import BottomBannerModalContent from "./ModalContent";
import { useConnect } from "./connect";
import { FORM_CONTAINER_ID } from "./constants";
import { Props } from "./types";

export default memo(function BottomBanner({
  includeSurveyLink = false,
  asFooter = false,
  customSurveyLinkPrefix,
  analyticsHandler,
  airtableLink,
  id = "newsletter-banner",
  isHubSpotReady: isHubSpotReadyProp = false,
  onHubSpotReady = noop,
}: Props): JSX.Element | null {
  const {
    setBottomBannerLastClosedTime,
    newsletterModalIsOpen,
    setIsHubSpotReady,
    isDirectLink,
    toggleNewsletterSignupModal,
    showBanner,
    isHubSpotReady,
    setEmail,
    setError,
    email,
    emailValidationError,
  } = useConnect({ isHubSpotReadyProp, asFooter });

  if (!showBanner) return null;

  return (
    <>
      <Head>
        {/* This is so that the mobile signup modal is scaled accordingly with the device. If using as footer there is not modal. */}
        {!asFooter && (
          <meta
            id="newsletter-signup-meta"
            name="viewport"
            content="width=device-width, initial-scale=1, minimum-scale=1, maximum-scale=5, user-scalable=yes"
          />
        )}
      </Head>
      {/* Script should be outside of <head> or else we get a next warning: "`next/script` should not be used in `next/head` component." */}
      <Script
        onReady={() => {
          setIsHubSpotReady(true);
          onHubSpotReady();
        }}
        type="text/javascript"
        src="https://js.hsforms.net/forms/v2.js"
      />
      <StyledBottomBannerWrapper
        asFooter={asFooter}
        id={BOTTOM_BANNER_ID}
        className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
        data-testid="newsletter-modal-banner-wrapper"
      >
        <StyledBanner
          dismissible={!asFooter}
          sdsType={"primary"}
          onClose={() => setBottomBannerLastClosedTime(Date.now())}
          data-id={id}
        >
          {/* Hidden form for submitting the data to HubSpot */}
          <HiddenHubSpotForm id={FORM_CONTAINER_ID} />

          {asFooter ? (
            <FooterContentWrapper>
              <BottomBannerModalContent
                isHubSpotReady={isHubSpotReady}
                setError={setError}
                setEmail={setEmail}
                email={email}
                emailValidationError={emailValidationError}
              />
            </FooterContentWrapper>
          ) : (
            <>
              <div>
                <StyledLink
                  onClick={toggleNewsletterSignupModal}
                  data-testid="newsletter-modal-open-button"
                >
                  Subscribe
                </StyledLink>{" "}
                to our newsletter to receive updates about new features.{" "}
                {includeSurveyLink && (
                  <>
                    {customSurveyLinkPrefix
                      ? customSurveyLinkPrefix
                      : "Send us feedback with this"}{" "}
                    <StyledLink
                      href={airtableLink}
                      target="_blank"
                      rel="noopener"
                      onClick={analyticsHandler}
                    >
                      quick survey
                    </StyledLink>
                    .
                  </>
                )}
              </div>

              <NewsletterModal
                isOpen={newsletterModalIsOpen}
                title=""
                onClose={toggleNewsletterSignupModal}
                isCloseButtonShown={false}
              >
                <HeaderContainer>
                  <Image alt="CellxGene Logo" src={cellxgeneLogoSvg} />
                  <StyledCloseButtonIcon
                    sdsIcon="xMark"
                    sdsSize="small"
                    onClick={toggleNewsletterSignupModal}
                    hideCloseButton={isDirectLink}
                    data-testid="newsletter-modal-close-button"
                  />
                </HeaderContainer>
                <BottomBannerModalContent
                  isHubSpotReady={isHubSpotReady}
                  setError={setError}
                  setEmail={setEmail}
                  email={email}
                  emailValidationError={emailValidationError}
                />
              </NewsletterModal>
            </>
          )}
        </StyledBanner>
      </StyledBottomBannerWrapper>
    </>
  );
});
