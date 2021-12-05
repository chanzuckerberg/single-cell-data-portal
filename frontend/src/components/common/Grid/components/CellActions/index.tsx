import React, { ReactNode } from "react";
import { Actions } from "src/components/common/Grid/components/CellActions/style";

interface Props {
  children: ReactNode;
}

export default function CellActions({ children }: Props): JSX.Element {
  return <Actions>{children}</Actions>;
}
