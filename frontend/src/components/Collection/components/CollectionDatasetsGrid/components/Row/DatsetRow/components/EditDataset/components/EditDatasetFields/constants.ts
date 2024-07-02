import { InputTextProps } from "@czi-sds/components";

type PickInputTextProps = "fullWidth" | "id" | "name" | "sdsType";

export const FIELD_NAMES = {
  TITLE: "title",
} as const;

export const INPUT_TEXT_PROPS: Record<
  keyof typeof FIELD_NAMES,
  Pick<InputTextProps, PickInputTextProps>
> = {
  TITLE: {
    fullWidth: true,
    id: "title",
    name: FIELD_NAMES.TITLE,
    sdsType: "textField",
  },
};
