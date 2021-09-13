import React, { FC } from "react";
import { RightAlignedDetailsCell } from "../../../common/style";
import { Skeleton } from "../common/Skeleton";

interface Props {
  cellCount: number | null;
  isLoading: boolean;
}

const CellCount: FC<Props> = ({ cellCount, isLoading }) => {
  if (isLoading) {
    return (
      <td>
        <Skeleton />
      </td>
    );
  }

  return <RightAlignedDetailsCell>{cellCount || "-"}</RightAlignedDetailsCell>;
};

export default CellCount;
