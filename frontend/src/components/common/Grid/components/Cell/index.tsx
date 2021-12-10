import React from "react";
import { CellProps } from "react-table";
import { DatasetRow } from "src/components/common/Filter/common/entities";

export default function Cell(
  props: CellProps<DatasetRow, string[]>
): JSX.Element[] {
  const {
    cell: { value },
  } = props;
  return value.map((v: string) => <div key={v}>{v}</div>);
}
