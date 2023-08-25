import Image from "next/image";
import Script from "next/script";
import { ReactElement, useEffect, useMemo, useRef, useState } from "react";
import { useLocalStorage } from "react-use";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "src/views/WheresMyGene/components/GeneSearchBar/components/SaveExport";
import {
  BOTTOM_BANNER_ID,
  NewsletterModal,
  StyledBanner,
  StyledBottomBannerWrapper,
  StyledDescription,
  StyledDisclaimer,
  StyledForm,
  StyledTitle,
  StyledInputText,
  StyledLink,
  StyledSubmitButton,
  HeaderContainer,
  StyledErrorMessage,
  HiddenHubspotForm,
  FooterContentWrapper,
  StyledCloseButtonIcon,
} from "./style";
import cellxgeneLogoSvg from "./CellxGene.svg";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import Head from "next/head";

export const FORM_CONTAINER_ID = "hubspot-form-container";
export const FORM_CONTAINER_ID_QUERY = `#${FORM_CONTAINER_ID}`;

export const SUBMIT_ISSUE_URL = "https://airtable.com/shrLwepDSEX1HI6bo";

export const FAILED_EMAIL_VALIDATION_STRING =
  "Please provide a valid email address.";

export interface Props {
  includeSurveyLink?: boolean;
  asFooter?: boolean;
  customSurveyLinkPrefix?: ReactElement;
  analyticsHandler?: () => void;
  airtableLink: string;
}
const BOTTOM_BANNER_EXPIRATION_TIME_MS = 30 * 24 * 60 * 60 * 1000; // 30 days
const BOTTOM_BANNER_LAST_CLOSED_TIME_KEY = "bottomBannerLastClosedTime";

export default function BottomBanner({
  includeSurveyLink = false,
  asFooter = false,
  customSurveyLinkPrefix,
  analyticsHandler,
  airtableLink,
}: Props): JSX.Element | null {
  const [bottomBannerLastClosedTime, setBottomBannerLastClosedTime] =
    useLocalStorage<number>(BOTTOM_BANNER_LAST_CLOSED_TIME_KEY, 0);

  const [newsletterModalIsOpen, setNewsletterModalIsOpen] = useState(false);
  const [isHubSpotReady, setIsHubSpotReady] = useState(false);
  const [isSubmitted, setIsSubmitted] = useState(false);
  const [email, setEmail] = useState("");
  const [emailValidationError, setError] = useState("");
  const [isDirectLink, setIsDirectLink] = useState(false);

  // For analytics if submit button was made enabled by user input
  const [submitButtonEnabledOnce, setSubmitButtonEnabledOnce] = useState(false);

  function toggleNewsletterSignupModal() {
    // Track when modal is opened
    if (!newsletterModalIsOpen) {
      track(EVENTS.NEWSLETTER_OPEN_MODAL_CLICKED);
    }

    setError("");
    setEmail("");
    setNewsletterModalIsOpen(!newsletterModalIsOpen);
  }

  // eslint-disable-next-line sonarjs/cognitive-complexity
  useEffect(() => {
    // Reads a query parameter from the URL to auto open the newsletter signup modal
    // Allows sharing of a URL to lead directly to the newsletter signup, specifically for conferences
    if (!asFooter && window) {
      const openModalParam = new URLSearchParams(window.location.search).get(
        "newsletter_signup"
      );

      if (openModalParam && openModalParam.toLowerCase() === "true") {
        setNewsletterModalIsOpen(true);
        setIsDirectLink(() => {
          if (!isDirectLink) {
            // Tracking this in setter function or else we get window not defined error
            track(EVENTS.NEWSLETTER_DIRECT_LINK_NAVIGATED);
          }
          return true;
        });
      }
    }

    // Observer to observe changes in the Hubspot embedded form, which is hidden from the user in order to use our own form view
    const observer = new MutationObserver((mutations) => {
      for (const mutation of mutations) {
        // Loop through all added nodes that were detected
        for (let i = 0; i < mutation.addedNodes.length; i++) {
          const node = mutation.addedNodes.item(i);

          // Submission success flow
          if (
            node?.textContent?.includes("Thank you for joining our newsletter.")
          ) {
            setIsSubmitted(true);
            setError("");

            track(EVENTS.NEWSLETTER_SIGNUP_SUCCESS);
          }

          // Hubspot email validation failure flow
          else if (
            node?.textContent?.includes("Please enter a valid email address.")
          ) {
            // HTML email validation may pass, but may not pass validation for Hubspot
            // ex. "ashintest_04252023_invalid_email@contractor.chanzuckerberg" does not validate with Hubspot but does with HTML email validation
            setError(FAILED_EMAIL_VALIDATION_STRING);

            track(EVENTS.NEWSLETTER_SIGNUP_FAILURE);
          }
        }
      }
    });

    if (isHubSpotReady) {
      hbspt.forms.create({
        region: "na1",
        portalId: "7272273",
        formId: "eb65b811-0451-414d-8304-7b9b6f468ce5",
        target: FORM_CONTAINER_ID_QUERY,
      });

      const form = document.querySelector(FORM_CONTAINER_ID_QUERY);
      if (form) {
        observer.observe(form, {
          childList: true,
          subtree: true,
        });
      }
    }

    return () => observer.disconnect();
  }, [asFooter, isDirectLink, isHubSpotReady]);

  const emailRef = useRef<HTMLInputElement | null>(null);

  // Validates if the email is valid or missing
  const validate = () => {
    const validityState = emailRef.current?.validity;
    if (validityState?.valueMissing || validityState?.typeMismatch) {
      setError(FAILED_EMAIL_VALIDATION_STRING);

      track(EVENTS.NEWSLETTER_SIGNUP_FAILURE);

      return false;
    }
    setError(""); // email validation passed, no error
    return true;
  };

  const handleSubmit = (event: React.FormEvent) => {
    track(EVENTS.NEWSLETTER_EMAIL_SUBMITTED);

    event.preventDefault();
    const isValid = validate();
    const form: HTMLFormElement | null = isValid
      ? document.querySelector(`${FORM_CONTAINER_ID_QUERY} form`)
      : null;

    if (!isValid || !form) {
      return;
    }

    const input = form.querySelector("input");
    if (!input) {
      return;
    }

    input.value = email;
    input.dispatchEvent(new Event("input", { bubbles: true }));

    form.submit();
  };

  const modalContent = (
    <div data-testid="newsletter-modal-content">
      <StyledTitle>Join Our Newsletter</StyledTitle>

      <StyledDescription>
        {isSubmitted
          ? "Thanks for subscribing!"
          : "Get a quarterly email with the latest CELLxGENE features and data."}
      </StyledDescription>

      <StyledForm onSubmit={handleSubmit} noValidate>
        {!isSubmitted && (
          <>
            <StyledInputText
              intent={emailValidationError ? "error" : "default"}
              inputRef={emailRef}
              placeholder={"Enter email address"}
              label={"Email"}
              hideLabel
              onChange={(event) => {
                if (emailValidationError) setError("");

                if (!submitButtonEnabledOnce) {
                  setSubmitButtonEnabledOnce(true);
                  track(EVENTS.NEWSLETTER_SUBSCRIBE_BUTTON_AVAILABLE);
                }

                setEmail(event.target.value);
              }}
              id={"email-input"}
              value={email}
              required
              type="email"
              inputProps={{ "data-testid": "newsletter-email-input" }}
            />
            <StyledSubmitButton
              type="submit"
              color="primary"
              name="subscribe"
              variant="contained"
              disableElevation
              disabled={!email}
              data-testid="newsletter-subscribe-button"
            >
              Subscribe
            </StyledSubmitButton>
          </>
        )}
      </StyledForm>

      <StyledErrorMessage data-testid="newsletter-validation-error-message">
        {emailValidationError}
      </StyledErrorMessage>

      <StyledDisclaimer>
        {isSubmitted
          ? 'To unsubscribe, click on the "Unsubscribe" link at the bottom of the newsletter.'
          : "Unsubscribe at any time."}
      </StyledDisclaimer>
    </div>
  );

  const showBanner = useMemo(() => {
    const show =
      !bottomBannerLastClosedTime ||
      Date.now() - bottomBannerLastClosedTime >
        BOTTOM_BANNER_EXPIRATION_TIME_MS;
    if (show && bottomBannerLastClosedTime) {
      setBottomBannerLastClosedTime(0);
    }
    return show;
  }, [bottomBannerLastClosedTime, setBottomBannerLastClosedTime]);

  if (!showBanner) return null;

  return (
    <>
      <Head>
        {/* This is so that the mobile signup modal is scaled accordingly with the device. If using as footer there is not modal. */}
        {!asFooter && (
          <meta
            id="newsletter-signup-meta"
            name="viewport"
            content="width=device-width, initial-scale=1, maximum-scale=1"
          />
        )}
      </Head>

      {/* Script should be outside of <head> or else we get a next warning: "`next/script` should not be used in `next/head` component." */}
      <Script
        onReady={() => {
          setIsHubSpotReady(true);
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
        >
          {/* Hidden form for submitting the data to Hubspot */}
          <HiddenHubspotForm id={FORM_CONTAINER_ID} />

          {asFooter ? (
            <FooterContentWrapper>{modalContent}</FooterContentWrapper>
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
                {modalContent}
              </NewsletterModal>
            </>
          )}
        </StyledBanner>
      </StyledBottomBannerWrapper>
    </>
  );
}
