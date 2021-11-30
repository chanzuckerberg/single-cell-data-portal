import loadable from "@loadable/component";
import { FC } from "react";
import { PluralizedMetadataLabel } from "src/common/constants/metadata";
import { LeftAlignedDetailsCell } from "../../../common/style";
import { Skeleton } from "../common/Skeleton";

const AsyncPopover = loadable(
  () =>
    /*webpackChunkName: 'CollectionRow/components/Popover' */ import(
      "src/components/Collections/components/Grid/components/Popover"
    )
);

interface Props {
  label: PluralizedMetadataLabel;
  values: string[];
  isLoading: boolean;
}

const Popover: FC<Props> = ({ label, values, isLoading }) => {
  if (isLoading) {
    return (
      <td>
        <Skeleton />
      </td>
    );
  }

  if (!values || values.length === 0) {
    return <LeftAlignedDetailsCell>-</LeftAlignedDetailsCell>;
  }

  return <AsyncPopover label={label} values={values} />;
};

export default Popover;
