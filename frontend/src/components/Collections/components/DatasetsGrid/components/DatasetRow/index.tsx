import loadable from "@loadable/component";
import React, { FC } from "react";
import { Dataset } from "src/common/entities";
import { useDatasetStatus } from "src/common/queries/datasets";
import { LeftAlignedDetailsCell } from "src/components/Collections/components/CollectionsGrid/components/CollectionRow/components/common/style";
import { aggregateDatasetsMetadata } from "src/components/Collections/components/CollectionsGrid/components/CollectionRow/components/common/utils";
import {
  RightAlignedDetailsCell,
  StyledRow,
} from "src/components/Collections/components/CollectionsGrid/components/CollectionRow/style";
import { StyledCell } from "src/components/Collections/components/CollectionsGrid/components/common/style";

interface Props {
  dataset: Dataset;
}

const AsyncPopover = loadable(
  () =>
    /*webpackChunkName: 'CollectionRow/components/Popover' */ import(
      "src/components/Collections/components/CollectionsGrid/components/CollectionRow/components/Popover"
    )
);

const conditionalPopover = (values: string[]) => {
  if (!values || values.length === 0) {
    return <LeftAlignedDetailsCell>-</LeftAlignedDetailsCell>;
  }

  return <AsyncPopover values={values} />;
};

const DatasetRow: FC<Props> = ({ dataset }) => {
  const {
    tissue,
    assay,
    disease,
    organism,
    cell_count,
  } = aggregateDatasetsMetadata([dataset]);

  const { name } = dataset;

  const datasetStatus = useDatasetStatus(dataset.id);

  return (
    <StyledRow>
      <StyledCell>{name}</StyledCell>
      {conditionalPopover(tissue)}
      {conditionalPopover(assay)}
      {conditionalPopover(disease)}
      {conditionalPopover(organism)}
      <RightAlignedDetailsCell>{cell_count || "-"}</RightAlignedDetailsCell>
      <RightAlignedDetailsCell>
        {/* {uploading ? uploadingStatus() : exploreButton()} */}
      </RightAlignedDetailsCell>
    </StyledRow>
  );
};

const uploadStatus = () => {
  return null;
};

const exploreButton = () => {
  return null;
};

export default DatasetRow;
