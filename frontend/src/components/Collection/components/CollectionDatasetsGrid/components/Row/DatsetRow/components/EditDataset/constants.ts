import { ButtonProps, DialogProps } from "@czi-sds/components";

export const CANCEL_BUTTON_PROPS: ButtonProps = {
  isAllCaps: false,
  sdsStyle: "minimal",
  sdsType: "secondary",
};

export const DIALOG_PROPS: Partial<DialogProps> = {
  PaperProps: { component: "form", noValidate: true },
  sdsSize: "s",
};

export const SAVE_BUTTON_PROPS: ButtonProps = {
  sdsStyle: "square",
  sdsType: "primary",
  type: "submit",
};
