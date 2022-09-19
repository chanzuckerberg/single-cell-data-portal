import React from "react";
import { TableCountSummary } from "src/components/common/Filter/common/entities";
import {
  CountAndTotal,
  HeaderCell as Cell,
} from "src/components/common/Grid/components/HeaderCell/style";

interface Props {
  label: string;
  tableCountSummary?: TableCountSummary;
}

export default function HeaderCell({
  label,
  tableCountSummary,
}: Props): JSX.Element {
  const rowCount = tableCountSummary?.row;
  const totalCount = tableCountSummary?.total;
  const countAndTotal =
    rowCount && totalCount ? `${rowCount} of ${totalCount}` : "";
  return (
    <Cell>
      <span>{label}</span>
      {countAndTotal && <CountAndTotal>{countAndTotal}</CountAndTotal>}
    </Cell>
  );
}
