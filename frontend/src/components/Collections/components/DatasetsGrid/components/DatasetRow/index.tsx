import { Button, Checkbox, Intent, Spinner } from "@blueprintjs/core";
import loadable from "@loadable/component";
import React, { FC } from "react";
import { Dataset, DatasetUploadStatus } from "src/common/entities";
import { useDatasetStatus } from "src/common/queries/datasets";
import { LeftAlignedDetailsCell } from "src/components/Collections/components/CollectionsGrid/components/CollectionRow/components/common/style";
import { aggregateDatasetsMetadata } from "src/components/Collections/components/CollectionsGrid/components/CollectionRow/components/common/utils";
import {
  RightAlignedDetailsCell,
  StyledRow,
} from "src/components/Collections/components/CollectionsGrid/components/CollectionRow/style";
import { StyledCell } from "src/components/Collections/components/CollectionsGrid/components/common/style";
import { UploadingFile } from "src/components/DropboxChooser";
import { UploadStatusContainer } from "./style";

interface Props {
  dataset: Dataset;
  checkHandler: (id: string) => void;
  file?: UploadingFile;
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

const DatasetRow: FC<Props> = ({ dataset, checkHandler, file }) => {
  const {
    tissue,
    assay,
    disease,
    organism,
    cell_count,
  } = aggregateDatasetsMetadata([dataset]);

  let { name } = dataset;
  if (!name) name = file?.name ?? dataset.id;

  const { data: datasetStatus, isError } = useDatasetStatus(dataset.id);
  if (isError) console.error(datasetStatus);
  if (!datasetStatus) return null;

  const isUploading = datasetStatus.upload_progress < 1;
  return (
    <StyledRow>
      <StyledCell>
        <Checkbox
          onChange={() => checkHandler(dataset.id)}
          style={{ marginBottom: "none" }}
        />
        {name}
      </StyledCell>
      {conditionalPopover(tissue)}
      {conditionalPopover(assay)}
      {conditionalPopover(disease)}
      {conditionalPopover(organism)}
      <RightAlignedDetailsCell>{cell_count || "-"}</RightAlignedDetailsCell>
      <RightAlignedDetailsCell>
        {isUploading ? (
          renderUploadStatus(datasetStatus)
        ) : (
          <Button intent={Intent.PRIMARY} outlined text="Explore" />
        )}
      </RightAlignedDetailsCell>
    </StyledRow>
  );
};

const renderUploadStatus = (datasetStatus: DatasetUploadStatus) => {
  return (
    <UploadStatusContainer>
      <Spinner
        intent={Intent.PRIMARY}
        value={datasetStatus.upload_progress}
        size={16}
      />

      {`${datasetStatus.upload_status} (${Math.round(
        datasetStatus.upload_progress * 100
      )}%)`}
    </UploadStatusContainer>
  );
};

export default DatasetRow;
