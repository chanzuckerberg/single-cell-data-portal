import { useEffect, useMemo, useRef, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  FAILED_EMAIL_VALIDATION_STRING,
  FORM_CONTAINER_ID_QUERY,
  HIDDEN_NEWSLETTER_SUCCESS_CLASS,
} from "../../constants";

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
  const emailRef = useRef<HTMLInputElement | null>(null);

  /**
   * handleSubmit
   * handles the submission of a newsletter signup form and includes
   * validation, verfication of elements, and event tracking
   */
  const handleSubmit = (event: React.FormEvent) => {
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

    track(EVENTS.NEWSLETTER_EMAIL_SUBMITTED);
  };

  /**
   * Validates the email input field by checking its validity state.
   *
   * - If the email field is missing a value or contains an invalid email format,
   *   an error message is set, and a failure event is tracked.
   * - If the email is valid, the error message is cleared.
   *
   * @returns {boolean} `true` if the email input is valid, otherwise `false`.
   */
  const validate = () => {
    const validityState = emailRef.current?.validity;

    if (validityState?.valueMissing || validityState?.typeMismatch) {
      setError(FAILED_EMAIL_VALIDATION_STRING);
      track(EVENTS.NEWSLETTER_SIGNUP_FAILURE);
      return false;
    }

    setError("");
    return true;
  };

  /**
   * formContainerQueryId
   * query selector for the form container
   */
  const formContainerQueryId = useMemo(() => {
    return `[data-id="${id}"] ${FORM_CONTAINER_ID_QUERY}`;
  }, [id]);

  /**
   * useEffect
   * handles the submission success flow and email validation failure flow
   */
  useEffect(() => {
    if (!isHubSpotReady) return;
    /**
     * Observer to observe changes in the Hubspot embedded form,
     * which is hidden from the user in order to use our own form view
     */
    const observer = new MutationObserver((mutations) => {
      for (const mutation of mutations) {
        /**
         * Loop through all added nodes that were detected
         */
        for (let i = 0; i < mutation.addedNodes.length; i++) {
          const node = mutation.addedNodes.item(i);
          /**
           *  Submission success flow
           */
          if (
            node instanceof Element &&
            node.classList.contains(HIDDEN_NEWSLETTER_SUCCESS_CLASS)
          ) {
            setIsSubmitted(true);
            setError("");
            track(EVENTS.NEWSLETTER_SIGNUP_SUCCESS);
          } else if (
            /**
             * Hubspot email validation failure flow
             */
            node?.textContent?.includes("Please enter a valid email address.")
          ) {
            /**
             * HTML email validation may pass, but may not pass validation for Hubspot
             */
            setError(FAILED_EMAIL_VALIDATION_STRING);
            track(EVENTS.NEWSLETTER_SIGNUP_FAILURE);
          }
        }
      }
    });

    const form = document.querySelector(formContainerQueryId);

    hbspt.forms.create({
      region: "na1",
      portalId: "7272273",
      formId: "eb65b811-0451-414d-8304-7b9b6f468ce5",
      target: formContainerQueryId,
      formInstanceId: id,
    });

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
