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
    formInstanceId?: string;
  }

  export interface HubSpotForm {
    create(options: CreateFormOptions): void;
  }
  export interface HubSpot {
    forms: HubSpotForm;
  }
}

declare const hbspt: HubSpotFormAPI.HubSpot;

declare module "src/components/common/Filter/descendant_mappings/cell_type_descendants.json" {
  const value: { [key: string]: string[] };
  export default value;
}

declare module "src/components/common/Filter/descendant_mappings/tissue_descendants.json" {
  const value: { [key: string]: string[] };
  export default value;
}
