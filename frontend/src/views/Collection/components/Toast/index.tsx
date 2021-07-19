import { IToaster, IToastProps, Position, Toaster } from "@blueprintjs/core";

let toaster: IToaster | null = null;

const Toast = {
  show(props: IToastProps) {
    if (!toaster) {
      toaster = Toaster.create({
        position: Position.TOP,
      });
    }

    toaster.show(props);
  },
};

export default Toast;
