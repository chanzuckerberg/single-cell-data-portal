import { API_URL } from "src/configs/configs";
import { EVENTS } from "./events";

export function track(event: EVENTS, props?: Record<string, unknown>): void {
  const options = props ? { props } : undefined;

  window.plausible(event, options);
  if (API_URL !== "https://api.cellxgene.cziscience.com")
    console.info(`Tracking event: ${event}`, options);
}
