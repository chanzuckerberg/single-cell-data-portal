import { EVENTS } from "./events";

export function track(event: EVENTS, props?: Record<string, unknown>): void {
  const options = props ? { props } : undefined;

  window.plausible(event, options);
}
