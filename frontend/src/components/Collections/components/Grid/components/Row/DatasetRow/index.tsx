import {
  Button,
  Checkbox,
  Classes,
  Icon,
  Intent,
  Spinner,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import loadable from "@loadable/component";
import React, { FC, useState } from "react";
import { useQueryCache } from "react-query";
import {
  Dataset,
  DatasetUploadStatus,
  UPLOAD_STATUS,
  VALIDATION_STATUS,
} from "src/common/entities";
import {
  useDatasetStatus,
  USE_DATASET_STATUS,
} from "src/common/queries/datasets";
import { aggregateDatasetsMetadata } from "src/components/Collections/components/Grid/common/utils";
import {
  DetailsCell,
  LeftAlignedDetailsCell,
  RightAlignedDetailsCell,
  StyledRow,
} from "src/components/Collections/components/Grid/components/Row/common/style";
import { UploadingFile } from "src/components/DropboxChooser";
import DatasetUploadToast from "src/views/Collection/components/DatasetUploadToast";
import { DatasetStatusTag, TitleContainer } from "./style";

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

const conditionalPopover = (values: string[], loading?: boolean) => {
  if (loading) return <td>{skeletonDiv}</td>;
  if (!values || values.length === 0) {
    return <LeftAlignedDetailsCell>-</LeftAlignedDetailsCell>;
  }

  return <AsyncPopover values={values} />;
};

const INITIAL_UPLOAD_PROGRESS = -1;

const DatasetRow: FC<Props> = ({ dataset, checkHandler, file }) => {
  const {
    tissue,
    assay,
    disease,
    organism,
    cell_count,
  } = aggregateDatasetsMetadata([dataset]);
  const queryCache = useQueryCache();

  let { name } = dataset;
  let datasetStatus = {} as DatasetUploadStatus;
  const queryResult = useDatasetStatus(dataset.id);
  const [lastUploadProgress, setLastUploadProgress] = useState(
    INITIAL_UPLOAD_PROGRESS
  );
  let statusFailed, isLoading;

  // TODO: When checking for conversion, will have to stop polling when conversion is done and there is no need for anymore checks

  // If there is no name on the dataset the conversion and upload process hasn't completed
  // Assign a temp name and begin polling the status endpoint
  // This should be replaced with a signifier from the backend instead of relying on name population
  if (!name) {
    name = file?.name ?? dataset.id;
    const { isError } = queryResult;
    if (isError) console.error(datasetStatus);
    if (!queryResult.data) return null;
    datasetStatus = queryResult.data;

    const isUploading =
      datasetStatus.upload_status === UPLOAD_STATUS.WAITING ||
      datasetStatus.upload_status === UPLOAD_STATUS.UPLOADING;
    const isPopulated = dataset.name !== "";
    statusFailed =
      datasetStatus.upload_status === UPLOAD_STATUS.FAILED ||
      datasetStatus.validation_status === VALIDATION_STATUS.INVALID;
    isLoading = !isPopulated || isUploading;

    if (statusFailed) {
      queryCache.cancelQueries([USE_DATASET_STATUS, dataset.id]);
      if (datasetStatus.upload_status === UPLOAD_STATUS.FAILED)
        DatasetUploadToast.show({
          intent: Intent.DANGER,
          message: "There was a problem uploading your file. Please try again.",
        });
    }
    if (lastUploadProgress !== datasetStatus.upload_progress) {
      if (
        datasetStatus.upload_progress === 1.0 &&
        lastUploadProgress !== INITIAL_UPLOAD_PROGRESS
      ) {
        DatasetUploadToast.show({
          icon: IconNames.TICK,
          intent: Intent.SUCCESS,
          message:
            "Upload was successful. Your file is being processed which will continue in the background, even if you close this window.",
        });
      }
      setLastUploadProgress(datasetStatus.upload_progress);
    }
  }
  return (
    <StyledRow>
      <DetailsCell>
        <TitleContainer>
          <Checkbox onChange={() => checkHandler(dataset.id)} />
          <div>{name}</div>
        </TitleContainer>
        {(isLoading || statusFailed) && renderUploadStatus(datasetStatus)}
      </DetailsCell>
      {conditionalPopover(tissue, isLoading)}
      {conditionalPopover(assay, isLoading)}
      {conditionalPopover(disease, isLoading)}
      {conditionalPopover(organism, isLoading)}
      {isLoading ? (
        <td>skeletonDiv</td>
      ) : (
        <RightAlignedDetailsCell>{cell_count}</RightAlignedDetailsCell>
      )}
      <RightAlignedDetailsCell>
        {!isLoading && (
          <Button intent={Intent.PRIMARY} outlined text="Explore" />
        )}
      </RightAlignedDetailsCell>
    </StyledRow>
  );
};

const renderUploadStatus = (datasetStatus: DatasetUploadStatus) => {
  if (datasetStatus.upload_status === UPLOAD_STATUS.FAILED)
    return (
      <DatasetStatusTag intent={Intent.DANGER}>
        <Icon iconSize={16} icon={IconNames.ISSUE} />
        Upload Error
      </DatasetStatusTag>
    );
  if (datasetStatus.upload_progress === 1.0)
    return (
      <DatasetStatusTag>
        <Spinner intent={Intent.PRIMARY} size={16} />
      </DatasetStatusTag>
    );

  return (
    <DatasetStatusTag>
      <Spinner
        intent={Intent.PRIMARY}
        value={datasetStatus.upload_progress}
        size={16}
      />

      {`${datasetStatus.upload_status} (${Math.round(
        datasetStatus.upload_progress * 100
      )}%)`}
    </DatasetStatusTag>
  );
};

export default DatasetRow;
