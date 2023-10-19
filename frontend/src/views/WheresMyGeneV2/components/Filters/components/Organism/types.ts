import {
  DefaultMenuSelectOption,
  InputDropdownProps as RawInputDropdownProps,
} from "@czi-sds/components";

export const InputDropdownProps: Partial<RawInputDropdownProps> = {
  sdsStyle: "square",
};

export interface Props {
  isLoading: boolean;
}

// (HACK): Not sure why styled Dropdown changes `onChange` type
export type tempOnChange = (
  options: DefaultMenuSelectOption | DefaultMenuSelectOption[] | null
) => void;
