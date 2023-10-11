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
  const emailRef = useRef<HTMLInputElement | null>(null);

  /**
   * handleSubmit
   * handles the submission of a newsletter signup form and includes
   * validation, verfication of elements, and event tracking
   */
  const handleSubmit = (event: React.FormEvent) => {
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
    event.preventDefault();
  };

  /**
   * validate
   * validates the email input field
   * returns true if the email is valid, false otherwise
   * sets the error state if the email is invalid
   * */
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

    const observer = new MutationObserver((mutations) => {
      for (const mutation of mutations) {
        for (let i = 0; i < mutation.addedNodes.length; i++) {
          const node = mutation.addedNodes.item(i);

          if (node?.textContent?.includes(HIDDEN_NEWSLETTER_SUCCESS_MESSAGE)) {
            setIsSubmitted(true);
            setError("");
            track(EVENTS.NEWSLETTER_SIGNUP_SUCCESS);
          } else if (
            node?.textContent?.includes("Please enter a valid email address.")
          ) {
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
