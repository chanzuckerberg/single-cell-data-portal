import React from "react";
import { TableCountSummary } from "src/components/common/Filter/common/entities";
import { CountAndTotal as Tag } from "./style";

interface Props {
  tableCountSummary?: TableCountSummary;
}

export default function CountAndTotal({
  tableCountSummary,
}: Props): JSX.Element {
  const rowCount = tableCountSummary?.row;
  const totalCount = tableCountSummary?.total;
  const countAndTotal =
    rowCount && totalCount ? `${rowCount} of ${totalCount}` : "";
  return (
    <>
      {!!countAndTotal && (
        <Tag
          color="info" // TOOD(SDSv20): This was gray
          hover={false}
          label={countAndTotal}
          sdsStyle="square"
          sdsType="primary"
        />
      )}
    </>
  );
}
