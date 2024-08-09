import { ReactElement } from "react";

export interface Props {
  hasSurveyLink?: boolean;
  hasNewsletterSignup?: boolean;
  asFooter?: boolean;
  customSurveyLinkPrefix?: ReactElement;
  analyticsHandler?: () => void;
  surveyLink: string;
  /**
   * (thuang): Since the homepage has two newsletter banners (footer and banner),
   * we need to give each banner a unique id to differentiate them.
   */
  id?: string;
  emailValidationError?: string;
}
