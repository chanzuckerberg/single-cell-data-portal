import React from "react";
import {
  CountAndTotal,
  HeaderCell as Cell,
} from "src/components/common/Grid/components/HeaderCell/style";

interface Props {
  label: string;
  rowCount?: number;
  totalCount?: number;
}

export default function HeaderCell({
  label,
  rowCount,
  totalCount,
}: Props): JSX.Element {
  const countAndTotal =
    rowCount && totalCount ? `${rowCount} of ${totalCount}` : "";
  return (
    <Cell>
      <span>{label}</span>
      {countAndTotal && <CountAndTotal>{countAndTotal}</CountAndTotal>}
    </Cell>
  );
}
