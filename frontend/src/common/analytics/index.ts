import { API_URL } from "src/configs/configs";
import { EVENTS } from "./events";

const serializeProps = (
  props?: Record<string, unknown>
): { props: Record<string, unknown> } | undefined => {
  if (!props) return undefined;
  return {
    props: Object.fromEntries(
      Object.entries(props).map(([key, val]) => {
        if (Array.isArray(val)) {
          const newSortedArray = [...val].sort();
          return [key as keyof typeof props, newSortedArray.join(",")];
        }
        return [key as keyof typeof props, val];
      })
    ),
  };
};

export function track(event: EVENTS, props?: Record<string, unknown>): void {
  const options = serializeProps(props);
  try {
    window.plausible(event, options);
  } catch (e) {
    // adding this to catch errors in tests
    console.error(e);
  }

  if (API_URL !== "https://api.cellxgene.cziscience.com") {
    console.info(`Tracking event: ${event}`, options?.props);
  }
}
