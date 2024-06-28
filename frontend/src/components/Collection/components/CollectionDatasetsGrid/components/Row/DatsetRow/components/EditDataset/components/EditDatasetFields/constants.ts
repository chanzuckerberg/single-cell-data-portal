import { InputTextProps } from "@czi-sds/components";

type PickInputTextProps = "fullWidth" | "id" | "name" | "sdsStage" | "sdsType";

export const TITLE_INPUT_TEXT_PROPS: Pick<InputTextProps, PickInputTextProps> =
  {
    fullWidth: true,
    id: "title",
    name: "title",
    sdsStage: "default",
    sdsType: "textField",
  };
