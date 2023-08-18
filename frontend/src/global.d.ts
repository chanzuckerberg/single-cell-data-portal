namespace HubSpotFormAPI {
  /**
   * HubSpot form creation options.
   * https://legacydocs.hubspot.com/docs/methods/forms/advanced_form_options
   */
  export interface CreateFormOptions {
    region: string;
    portalId: string;
    formId: string;
    target?: string;
  }

  export interface HubSpotForm {
    create(options: CreateFormOptions): void;
  }
  export interface HubSpot {
    forms: HubSpotForm;
  }
}

declare const hbspt: HubSpotFormAPI.HubSpot;
