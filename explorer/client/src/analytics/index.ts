import { EVENTS } from "./events";

declare global {
  interface Window {
    plausible: {
      q: unknown[];
      (event: EVENTS, options?: { props: { [key: string]: unknown } }): void;
    };
  }
}

export function track(event: EVENTS, props?: Record<string, unknown>): void {
  const options = props ? { props } : undefined;

  try {
    window.plausible(event, options);
  } catch (error) {
    console.error(error);
  }
}
