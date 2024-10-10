import { useEffect, useMemo, useState } from "react";
import { useLocalStorage } from "react-use";
import {
  BOTTOM_BANNER_EXPIRATION_TIME_MS,
  BOTTOM_BANNER_LAST_CLOSED_TIME_KEY,
} from "./constants";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";

export const useConnect = ({
  isHubSpotReadyProp,
  asFooter,
}: {
  isHubSpotReadyProp: boolean;
  asFooter?: boolean;
}) => {
  const [bottomBannerLastClosedTime, setBottomBannerLastClosedTime] =
    useLocalStorage<number>(BOTTOM_BANNER_LAST_CLOSED_TIME_KEY, 0);

  const [newsletterModalIsOpen, setNewsletterModalIsOpen] = useState(false);
  const [isHubSpotReady, setIsHubSpotReady] = useState(false);
  const [email, setEmail] = useState("");
  const [emailValidationError, setError] = useState("");
  const [isDirectLink, setIsDirectLink] = useState(false);

  useEffect(() => {
    if (!isHubSpotReadyProp) return;

    setIsHubSpotReady(true);
  }, [isHubSpotReady, isHubSpotReadyProp]);

  /**
   * useEffect
   * Reads a query parameter from the URL to auto open the newsletter signup modal
   * Allows sharing of a URL to lead directly to the newsletter signup, specifically
   * for conferences
   * isDirectLink is tracked in setter function or else we get window not defined error
   */
  useEffect(() => {
    if (!isHubSpotReady) return;

    if (!asFooter && window) {
      const openModalParam = new URLSearchParams(window.location.search).get(
        "newsletter_signup"
      );

      if (openModalParam && openModalParam.toLowerCase() === "true") {
        setNewsletterModalIsOpen(true);
        setIsDirectLink(() => {
          if (!isDirectLink) {
            track(EVENTS.NEWSLETTER_DIRECT_LINK_NAVIGATED);
          }
          return true;
        });
      }
    }
  }, [asFooter, isDirectLink, isHubSpotReady]);

  /**
   * showBanner
   * Returns true if the banner should be shown
   * If the banner has been closed in the last BOTTOM_BANNER_EXPIRATION_TIME_MS days, it
   * will not be shown. If the banner has never been closed, it will be shown
   */
  const showBanner = useMemo(() => {
    if (asFooter) return true;

    const show =
      !bottomBannerLastClosedTime ||
      Date.now() - bottomBannerLastClosedTime >
        BOTTOM_BANNER_EXPIRATION_TIME_MS;
    if (show && bottomBannerLastClosedTime) {
      setBottomBannerLastClosedTime(0);
    }
    return show;
  }, [asFooter, bottomBannerLastClosedTime, setBottomBannerLastClosedTime]);

  /**
   * toggleNewsletterSignupModal
   * Toggles the newsletter signup modal
   * */
  function toggleNewsletterSignupModal() {
    if (!newsletterModalIsOpen) {
      track(EVENTS.NEWSLETTER_OPEN_MODAL_CLICKED);
    }
    setError("");
    setEmail("");
    setNewsletterModalIsOpen(!newsletterModalIsOpen);
  }

  return {
    setBottomBannerLastClosedTime,
    newsletterModalIsOpen,
    setIsHubSpotReady,
    isDirectLink,
    toggleNewsletterSignupModal,
    showBanner,
    isHubSpotReady,
    email,
    setEmail,
    emailValidationError,
    setError,
  };
};
