import { InputDropdownProps as IInputDropdownProps } from "@czi-sds/components";

export interface Props {
  setIsScaled: (prevIsScaled: boolean) => void;
}

export const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = {
  sdsStyle: "square",
};
