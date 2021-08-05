import { IToaster, IToastProps, Position, Toaster } from "@blueprintjs/core";

/**
 * (thuang): We can't lazy-load toaster via `.create()` on the first `.show()`
 * call, because sometimes `.show()` is called inside React lifecycle hooks.
 */
let toaster: IToaster | null = null;

// (thuang): Instantiates on the client side only
if (typeof document !== "undefined") {
  toaster = Toaster.create({
    position: Position.TOP,
  });
}

const Toast = {
  show(props: IToastProps): void {
    toaster?.show(props);
  },
};

export default Toast;
