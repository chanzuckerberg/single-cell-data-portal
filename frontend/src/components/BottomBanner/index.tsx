import { memo } from "react";
import {
  BOTTOM_BANNER_ID,
  NewsletterModal,
  StyledBanner,
  StyledBottomBannerWrapper,
  StyledLink,
  FooterContentWrapper,
  StyledCloseButtonIcon,
  HeaderContainer,
} from "./style";
import CellxgeneLogoSvg from "src/common/images/CellxGene.svg";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SaveExport";
import BottomBannerModalContent from "./components/ModalContent";
import { useConnect } from "./connect";
import {
  BOTTOM_BANNER_SURVEY_LINK_TEXT,
  BOTTOM_BANNER_SURVEY_TEXT,
  NEWSLETTER_SIGNUP_BANNER_SUBSCRIBE_BUTTON_TEXT,
  NEWSLETTER_SIGNUP_BANNER_SUBSCRIBE_TEXT,
} from "./constants";
import { Props } from "./types";

export default memo(function BottomBanner({
  hasSurveyLink = true,
  hasNewsletterSignup = false,
  asFooter = false,
  customSurveyLinkPrefix,
  analyticsHandler,
  surveyLink,
  id = "newsletter-banner",
}: Props): JSX.Element | null {
  const {
    setBottomBannerLastClosedTime,
    setEmail,
    setError,
    toggleNewsletterSignupModal,
    newsletterModalIsOpen,
    isDirectLink,
    showBanner,
    email,
    emailValidationError,
  } = useConnect({ asFooter });

  if (!showBanner) return null;

  return (
    <>
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
          {asFooter ? (
            <FooterContentWrapper>
              <BottomBannerModalContent
                setError={setError}
                setEmail={setEmail}
                email={email}
                emailValidationError={emailValidationError}
              />
            </FooterContentWrapper>
          ) : (
            <>
              <>
                {hasNewsletterSignup && (
                  <>
                    <StyledLink
                      onClick={toggleNewsletterSignupModal}
                      data-testid="newsletter-modal-open-button"
                    >
                      {NEWSLETTER_SIGNUP_BANNER_SUBSCRIBE_BUTTON_TEXT}
                    </StyledLink>
                    {NEWSLETTER_SIGNUP_BANNER_SUBSCRIBE_TEXT}
                  </>
                )}
                {hasSurveyLink && (
                  <>
                    {customSurveyLinkPrefix
                      ? customSurveyLinkPrefix
                      : BOTTOM_BANNER_SURVEY_TEXT}
                    <StyledLink
                      href={surveyLink}
                      target="_blank"
                      rel="noopener"
                      onClick={analyticsHandler}
                    >
                      {BOTTOM_BANNER_SURVEY_LINK_TEXT}
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
                  <CellxgeneLogoSvg />
                  <StyledCloseButtonIcon
                    icon="XMark"
                    sdsSize="small"
                    onClick={toggleNewsletterSignupModal}
                    hideCloseButton={isDirectLink}
                    data-testid="newsletter-modal-close-button"
                  />
                </HeaderContainer>
                <BottomBannerModalContent
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
