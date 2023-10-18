import { ReactElement } from "react";

export interface Props {
  includeSurveyLink?: boolean;
  asFooter?: boolean;
  customSurveyLinkPrefix?: ReactElement;
  analyticsHandler?: () => void;
  airtableLink: string;
  /**
   * (thuang): Since the homepage has two newsletter banners (footer and banner),
   * we need to give each banner a unique id to differentiate them.
   */
  id?: string;
  /**
   * (thuang): This is needed for homepage, because we have two instances of BottomBanner
   * and the HubSpot script is only loaded once, so only the first instance of BottomBanner
   * gets the onReady callback.
   */
  isHubSpotReady?: boolean;
  /**
   * (thuang): The first instance of BottomBanner in the UI needs to call onHubSpotReady(),
   * so the second BottomBanner knows when HubSpot is ready via `isHubSpotReady` prop.
   */
  onHubSpotReady?: () => void;
  emailValidationError?: string;
}
