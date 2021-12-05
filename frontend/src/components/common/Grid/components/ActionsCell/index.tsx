import React, { ReactNode } from "react";
import { Actions } from "src/components/common/Grid/components/ActionsCell/style";

interface Props {
  children: ReactNode;
}

export default function ActionsCell({ children }: Props): JSX.Element {
  return <Actions>{children}</Actions>;
}
