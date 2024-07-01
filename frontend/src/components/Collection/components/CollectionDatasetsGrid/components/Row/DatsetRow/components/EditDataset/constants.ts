import { ButtonProps, DialogProps } from "@czi-sds/components";

export const CANCEL_BUTTON_PROPS: ButtonProps = {
  isAllCaps: false,
  sdsStyle: "minimal",
  sdsType: "secondary",
};

export const DATASET_EDIT_CANCEL = "dataset-edit-cancel";
export const DATASET_EDIT_FORM = "dataset-edit-form";
export const DATASET_EDIT_SAVE = "dataset-edit-save";

export const DIALOG_PROPS: Partial<DialogProps> = {
  PaperProps: {
    component: "form",
    "data-testid": DATASET_EDIT_FORM,
    noValidate: true,
  },
  sdsSize: "s",
};

export const SAVE_BUTTON_PROPS: ButtonProps = {
  sdsStyle: "square",
  sdsType: "primary",
  type: "submit",
};
