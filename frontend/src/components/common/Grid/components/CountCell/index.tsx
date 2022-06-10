import React from "react";
import { CellCount } from "src/components/common/Grid/components/CountCell/style";

interface Props {
  cellCount: number;
}

export default function CountCell({ cellCount }: Props): JSX.Element {
  const count = cellCount ? cellCount.toLocaleString() : "";
  return <CellCount>{count}</CellCount>;
}
