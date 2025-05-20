import { ReactNode } from "react";
import { Caption } from "./style";

export default function CopyCaption({
  children,
}: {
  children: ReactNode;
}): JSX.Element {
  return <Caption>{children}</Caption>;
}
