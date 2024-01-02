import { createRoot } from "react-dom/client";

import {
  ToasterInstance,
  ToastProps,
  Position,
  OverlayToaster,
} from "@blueprintjs/core";

/**
 * (thuang): We can't lazy-load toaster on the first `.show()` call, because
 * sometimes `.show()` is called inside React lifecycle hooks.
 */
let toaster: ToasterInstance | null = null;

/**
 * (thuang): Instantiates on the client side only
 */
if (typeof document !== "undefined") {
  /**
   * (thuang): Don't use OverlayToaster.create() because it will use React 17
   * render() API, which is not supported by React 18
   * BluePrint plans to fix this in v6
   * https://github.com/palantir/blueprint/issues/5212#issuecomment-1124544078
   *
   * NOTE: In dev mode, when updating this file you might see warnings like:
   * "Warning: You are calling ReactDOMClient.createRoot() on a container that
   * has already been passed to createRoot() before. Instead, call root.render()
   * on the existing root instead if you want to update it."
   * This is due to hot module replacement (HMR) in dev mode. It's safe to ignore
   */
  createRoot(
    /**
     * (thuang): <div id="bp-toaster" /> is rendered in _app.tsx
     */
    // eslint-disable-next-line @typescript-eslint/no-non-null-assertion
    document.getElementById("bp-toaster")!
  ).render(
    <OverlayToaster
      position={Position.TOP}
      ref={(instance) => {
        toaster = instance;
      }}
    />
  );
}

const Toast = {
  show(props: ToastProps): void {
    toaster?.show(props);
  },
};

export default Toast;
