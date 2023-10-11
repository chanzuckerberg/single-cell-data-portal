import { AnyAction } from "redux";
import type { Config } from "../globals";

interface ConfigState {
  features: Config["features"] | null;
  parameters: Config["parameters"] | null;
  displayNames: Config["displayNames"] | null;
  loading?: boolean;
  error?: Error | null;
  corpora_props?: Config["corpora_props"];
  library_versions?: Config["library_versions"];
  portalUrl?: Config["portalUrl"];
  links?: Config["links"];
}
const ConfigReducer = (
  state: ConfigState = {
    displayNames: null,
    features: null,
    parameters: null,
  },
  action: AnyAction,
) => {
  switch (action.type) {
    case "initial data load start":
      return {
        ...state,
        loading: true,
        error: null,
      };
    case "configuration load complete":
      return {
        ...state,
        loading: false,
        error: null,
        ...action.config,
      };
    case "initial data load error":
      return {
        ...state,
        error: action.error,
      };
    default:
      return state;
  }
};

export default ConfigReducer;
