import { ButtonIcon } from "czifui";
import Image from "next/image";
import Script from "next/script";
import { useEffect, useRef, useState } from "react";
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
} from "./style";
import cellxgeneLogoSvg from "./CellxGene.svg";

export const FORM_CONTAINER_ID = "hubspot-form-container";
export const FORM_CONTAINER_ID_QUERY = `#${FORM_CONTAINER_ID}`;

export const SUBMIT_ISSUE_URL = "https://airtable.com/shrLwepDSEX1HI6bo";

export interface Props {
  survey?: boolean;
  asFooter?: boolean;
}

export default function BottomBanner({
  survey = false,
  asFooter = false,
}: Props): JSX.Element {
  const [newsletterModalIsOpen, setNewsletterModalIsOpen] = useState(false);
  const [isHubSpotReady, setIsHubSpotReady] = useState(false);
  const [isSubmitted, setIsSubmitted] = useState(false);
  const [email, setEmail] = useState("");
  const [error, setError] = useState("");

  function toggleNewsletterSignupModal() {
    setError("");
    setEmail("");
    setNewsletterModalIsOpen(!newsletterModalIsOpen);
  }

  useEffect(() => {
    // Reads a query parameter from the URL to auto open the newsletter signup modal
    // Allows sharing of a URL to lead directly to the newsletter signup, specifically for conferences
    if (window) {
      const openModalParam = new URLSearchParams(window.location.search).get(
        "newsletter_signup"
      );

      if (openModalParam) {
        setNewsletterModalIsOpen(openModalParam.toLowerCase() === "true");
      }
    }

    // Observer to observe changes in the Hubspot embedded form, which is hidden from the user in order to use our own form view
    const observer = new MutationObserver((mutations) => {
      for (const mutation of mutations) {
        for (let i = 0; i < mutation.addedNodes.length; i++) {
          const node = mutation.addedNodes.item(i);
          if (
            node?.textContent?.includes("Thank you for joining our newsletter.")
          ) {
            console.log("is submitted!!!!");
            setIsSubmitted(true);
          }
        }
      }
    });

    if (isHubSpotReady) {
      console.log("Creating Hubspot form. Target: " + FORM_CONTAINER_ID_QUERY);
      hbspt.forms.create({
        region: "na1",
        portalId: "7272273",
        formId: "eb65b811-0451-414d-8304-7b9b6f468ce5",
        target: FORM_CONTAINER_ID_QUERY,
      });

      const form = document.querySelector(FORM_CONTAINER_ID_QUERY);
      if (form) {
        console.log("observing...");
        observer.observe(form, {
          childList: true,
          subtree: true,
        });
      }
    }

    return () => observer.disconnect();
  }, [isHubSpotReady]);

  const emailRef = useRef<HTMLInputElement | null>(null);

  const validate = () => {
    const validityState = emailRef.current?.validity;
    if (validityState?.valueMissing) {
      console.log("missingEmail");
      setError("missingEmail");
      return false;
    }
    if (validityState?.typeMismatch) {
      console.log("invalidEmail");
      setError("invalidEmail");
      return false;
    }
    console.log("no email error");
    setError(""); // no error
    return true;
  };

  const handleSubmit = (event: React.FormEvent) => {
    event.preventDefault();
    const isValid = validate();
    const form: HTMLFormElement | null = isValid
      ? document.querySelector(`${FORM_CONTAINER_ID_QUERY} form`)
      : null;

    if (!isValid || !form) {
      console.log("form not valid!!!");
      return;
    }

    const input = form.querySelector("input");
    if (!input) {
      console.log("input not valid!!!");
      return;
    }

    input.value = email;
    input.dispatchEvent(new Event("input", { bubbles: true }));

    form.submit();
  };

  // ashintest_04192023@contractor.chanzuckerberg.com

  const modalContent = (
    <>
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
              intent={error ? "error" : "default"}
              inputRef={emailRef}
              placeholder={"Enter email address"}
              label={"Email"}
              hideLabel
              onChange={(event) => {
                if (error) setError("");
                setEmail(event.target.value);
              }}
              id={"email-input"}
              value={email}
              required
              type="email"
            />
            <StyledSubmitButton
              type="submit"
              color="primary"
              name="subscribe"
              variant="contained"
              disableElevation
              data-testid="newsletter-submit-button"
              disabled={!email}
            >
              Subscribe
            </StyledSubmitButton>
          </>
        )}
      </StyledForm>

      <StyledErrorMessage>
        {error ? "Please provide a valid email address." : ""}
      </StyledErrorMessage>

      <StyledDisclaimer>
        {isSubmitted
          ? 'To unsubscribe, click on the "Unsubscribe" link at the bottom of the newsletter.'
          : "Unsubscribe at any time."}
      </StyledDisclaimer>
    </>
  );

  return (
    <>
      <Script
        onReady={() => {
          console.log("Script loaded");
          setIsHubSpotReady(true);
        }}
        type="text/javascript"
        src="https://js.hsforms.net/forms/v2.js"
      />

      <StyledBottomBannerWrapper
        asFooter={asFooter}
        id={BOTTOM_BANNER_ID}
        className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
      >
        <StyledBanner dismissible={!asFooter} sdsType={"primary"}>
          {/* Hidden form for submitting the data to Hubspot */}
          <HiddenHubspotForm id={FORM_CONTAINER_ID} />

          {asFooter ? (
            <FooterContentWrapper>{modalContent}</FooterContentWrapper>
          ) : (
            <div>
              <StyledLink onClick={toggleNewsletterSignupModal}>
                Subscribe
              </StyledLink>{" "}
              to our newsletter to receive updates about new features.{" "}
              {survey && (
                <>
                  Send us feedback with this{" "}
                  <StyledLink
                    href="https://airtable.com/shrLwepDSEX1HI6bo"
                    target="_blank"
                    rel="noopener"
                  >
                    quick survey
                  </StyledLink>
                  .
                </>
              )}
            </div>
          )}
        </StyledBanner>
      </StyledBottomBannerWrapper>

      <NewsletterModal
        isOpen={newsletterModalIsOpen}
        title=""
        onClose={toggleNewsletterSignupModal}
        isCloseButtonShown={false}
      >
        <HeaderContainer>
          <Image alt="CellxGene Logo" src={cellxgeneLogoSvg} />
          <ButtonIcon
            sdsIcon="xMark"
            sdsSize="small"
            onClick={toggleNewsletterSignupModal}
          />
        </HeaderContainer>
        {modalContent}
      </NewsletterModal>
    </>
  );
}
