import { IconProps, MenuProps as SDSProps } from "@czi-sds/components";

export const MENU_PROPS: Partial<SDSProps> = {
  anchorOrigin: {
    horizontal: "right",
    vertical: "bottom",
  },
  transformOrigin: { horizontal: "right", vertical: "top" },
};

export const DELETE_ICON_PROPS: IconProps<"TrashCan"> = {
  color: "red",
  sdsIcon: "TrashCan",
  sdsSize: "xs",
  sdsType: "static",
  shade: 400,
};

export const EDIT_ICON_PROPS: IconProps<"Edit"> = {
  color: "gray",
  sdsIcon: "Edit",
  sdsSize: "xs",
  sdsType: "static",
  shade: 400,
};
