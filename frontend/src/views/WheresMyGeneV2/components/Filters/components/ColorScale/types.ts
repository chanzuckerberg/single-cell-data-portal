import { InputDropdownProps as IInputDropdownProps } from "@czi-sds/components";
import { SORT_BY } from "src/views/WheresMyGeneV2/common/types";

export interface Props {
  setIsScaled: (prevIsScaled: boolean) => void;
}

export const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = {
  sdsStyle: "square",
};

export type ColorScaleOptionType = { id?: SORT_BY; name: string };
