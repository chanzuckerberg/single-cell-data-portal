import React from "react";
import { CellPropsValue } from "src/components/common/Filter/common/entities";

export default function CountCell(value: CellPropsValue): JSX.Element {
  const count = value.value?.toLocaleString() || "";
  return <div>{count}</div>;
}
