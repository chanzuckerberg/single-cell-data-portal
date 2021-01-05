import { Button, Checkbox, Classes, Intent, Spinner } from "@blueprintjs/core";
import loadable from "@loadable/component";
import React, { FC } from "react";
import { Dataset, DatasetUploadStatus } from "src/common/entities";
import { useDatasetStatus } from "src/common/queries/datasets";
import { aggregateDatasetsMetadata } from "src/components/Collections/components/Grid/common/utils";
import {
  DetailsCell,
  LeftAlignedDetailsCell,
  RightAlignedDetailsCell,
  StyledRow,
} from "src/components/Collections/components/Grid/components/Row/common/style";
import { UploadingFile } from "src/components/DropboxChooser";
import { TitleContainer, UploadStatusContainer } from "./style";

interface Props {
  dataset: Dataset;
  checkHandler: (id: string) => void;
  file?: UploadingFile;
}

const AsyncPopover = loadable(
  () =>
    /*webpackChunkName: 'CollectionRow/components/Popover' */ import(
      "src/components/Collections/components/Grid/components/Popover"
    )
);

const skeletonDiv = <div className={Classes.SKELETON}>PLACEHOLDER_TEXT</div>;

const conditionalPopover = (values: string[], loading: boolean) => {
  if (loading) return <td>{skeletonDiv}</td>;
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
  const isPopulated = dataset.name !== "";
  const isLoading = !isPopulated || isUploading;
  return (
    <StyledRow>
      <DetailsCell>
        <TitleContainer>
          <Checkbox onChange={() => checkHandler(dataset.id)} />
          <div>{name}</div>
        </TitleContainer>
      </DetailsCell>
      {conditionalPopover(tissue, isLoading)}
      {conditionalPopover(assay, isLoading)}
      {conditionalPopover(disease, isLoading)}
      {conditionalPopover(organism, isLoading)}
      <RightAlignedDetailsCell>
        {isLoading ? skeletonDiv : cell_count}
      </RightAlignedDetailsCell>
      <RightAlignedDetailsCell>
        {isLoading ? (
          renderUploadStatus(datasetStatus)
        ) : (
          <Button intent={Intent.PRIMARY} outlined text="Explore" />
        )}
      </RightAlignedDetailsCell>
    </StyledRow>
  );
};

const renderUploadStatus = (datasetStatus: DatasetUploadStatus) => {
  if (datasetStatus.upload_progress === 1.0)
    return (
      <UploadStatusContainer>
        <Spinner intent={Intent.PRIMARY} size={16} />
      </UploadStatusContainer>
    );

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
