import { track } from "src/common/analytics";
import {
  StyledDescription,
  StyledTitle,
  StyledForm,
  StyledInputText,
  StyledSubmitButton,
  StyledErrorMessage,
  StyledDisclaimer,
  StyledNewsletterSignupModal,
} from "./style";
import { useConnect } from "./connect";
import { EVENTS } from "src/common/analytics/events";
import { Props } from "./types";

export default function BottomBannerModalContent({
  id = "newsletter-banner",
  isHubSpotReady = false,
  setEmail,
  setError,
  emailValidationError,
  email,
}: Props): JSX.Element {
  const {
    isSubmitted,
    handleSubmit,
    submitButtonEnabledOnce,
    setSubmitButtonEnabledOnce,
    emailRef,
  } = useConnect({ id, isHubSpotReady, setError, email });
  return (
    <StyledNewsletterSignupModal data-testid="newsletter-modal-content">
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
    </StyledNewsletterSignupModal>
  );
}
