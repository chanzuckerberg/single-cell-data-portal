import { useEffect, useMemo, useRef, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  FAILED_EMAIL_VALIDATION_STRING,
  FORM_CONTAINER_ID_QUERY,
  HIDDEN_NEWSLETTER_SUCCESS_MESSAGE,
} from "../constants";

export const useConnect = ({
  id,
  isHubSpotReady,
  email,
  setError,
}: {
  id?: string;
  isHubSpotReady: boolean;
  email: string;
  setError: (error: string) => void;
}) => {
  const [isSubmitted, setIsSubmitted] = useState(false);
  const [submitButtonEnabledOnce, setSubmitButtonEnabledOnce] = useState(false);

  const handleSubmit = (event: React.FormEvent) => {
    track(EVENTS.NEWSLETTER_EMAIL_SUBMITTED);

    event.preventDefault();
    const isValid = validate();
    const form: HTMLFormElement | null = isValid
      ? document.querySelector(`${formContainerQueryId} form`)
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
  const formContainerQueryId = useMemo(() => {
    return `[data-id="${id}"] ${FORM_CONTAINER_ID_QUERY}`;
  }, [id]);

  useEffect(() => {
    if (!isHubSpotReady) return;

    // Observer to observe changes in the HubSpot embedded form, which is hidden from the user in order to use our own form view
    const observer = new MutationObserver((mutations) => {
      for (const mutation of mutations) {
        // Loop through all added nodes that were detected
        for (let i = 0; i < mutation.addedNodes.length; i++) {
          const node = mutation.addedNodes.item(i);

          // Submission success flow
          if (node?.textContent?.includes(HIDDEN_NEWSLETTER_SUCCESS_MESSAGE)) {
            setIsSubmitted(true);
            setError("");

            track(EVENTS.NEWSLETTER_SIGNUP_SUCCESS);
          }

          // HubSpot email validation failure flow
          else if (
            node?.textContent?.includes("Please enter a valid email address.")
          ) {
            // HTML email validation may pass, but may not pass validation for HubSpot
            // ex. "ashintest_04252023_invalid_email@contractor.chanzuckerberg" does not validate with HubSpot but does with HTML email validation
            setError(FAILED_EMAIL_VALIDATION_STRING);

            track(EVENTS.NEWSLETTER_SIGNUP_FAILURE);
          }
        }
      }
    });

    hbspt.forms.create({
      region: "na1",
      portalId: "7272273",
      formId: "eb65b811-0451-414d-8304-7b9b6f468ce5",
      target: formContainerQueryId,
      /**
       * (thuang): Since the homepage has two newsletter banners (footer and banner),
       * we need to give each banner a unique id to differentiate them.
       * https://legacydocs.hubspot.com/docs/methods/forms/advanced_form_options
       */
      formInstanceId: id,
    });

    const form = document.querySelector(formContainerQueryId);

    if (form) {
      observer.observe(form, {
        childList: true,
        subtree: true,
      });
    }

    return () => {
      observer.disconnect();
    };
  }, [formContainerQueryId, id, setError, isHubSpotReady]);

  return {
    isSubmitted,
    handleSubmit,
    emailRef,
    submitButtonEnabledOnce,
    setSubmitButtonEnabledOnce,
  };
};
