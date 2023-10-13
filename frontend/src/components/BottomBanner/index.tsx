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
import {
  FORM_CONTAINER_ID,
  HUBSPOT_URL,
  NEWSLETTER_SIGNUP_BANNER_SUBSCRIBE_BUTTON_TEXT,
  NEWSLETTER_SIGNUP_BANNER_SUBSCRIBE_TEXT,
  NEWSLETTER_SIGNUP_BANNER_SURVEY_LINK_TEXT,
  NEWSLETTER_SIGNUP_BANNER_SURVEY_TEXT,
} from "./constants";
import { Props } from "./types";

export default memo(function BottomBanner({
  includeSurveyLink = true,
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
    setEmail,
    setError,
    setIsHubSpotReady,
    toggleNewsletterSignupModal,
    newsletterModalIsOpen,
    isDirectLink,
    showBanner,
    isHubSpotReady,
    email,
    emailValidationError,
  } = useConnect({ isHubSpotReadyProp, asFooter });

  if (!showBanner) return null;

  return (
    <>
      <Head>
        {!asFooter && (
          <meta
            id="newsletter-signup-meta"
            name="viewport"
            content="width=device-width, initial-scale=1, maximum-scale=1"
          />
        )}
      </Head>
      <Script
        onReady={() => {
          setIsHubSpotReady(true);
          onHubSpotReady();
        }}
        type="text/javascript"
        src={HUBSPOT_URL}
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
              <>
                <StyledLink
                  onClick={toggleNewsletterSignupModal}
                  data-testid="newsletter-modal-open-button"
                >
                  {NEWSLETTER_SIGNUP_BANNER_SUBSCRIBE_BUTTON_TEXT}
                </StyledLink>
                {NEWSLETTER_SIGNUP_BANNER_SUBSCRIBE_TEXT}
                {includeSurveyLink && (
                  <>
                    {customSurveyLinkPrefix
                      ? customSurveyLinkPrefix
                      : NEWSLETTER_SIGNUP_BANNER_SURVEY_TEXT}
                    <StyledLink
                      href={airtableLink}
                      target="_blank"
                      rel="noopener"
                      onClick={analyticsHandler}
                    >
                      {NEWSLETTER_SIGNUP_BANNER_SURVEY_LINK_TEXT}
                    </StyledLink>
                  </>
                )}
              </>

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
