export interface Props {
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
  setError: (error: string) => void;
  setEmail: (email: string) => void;
  email: string;
  emailValidationError?: string;
}
