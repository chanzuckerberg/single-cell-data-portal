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
import { QueryCache, useQueryCache } from "react-query";
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

const conditionalPopover = (values: string[], loading?: boolean, hasFailed) => {
  if (loading) return <td>{skeletonDiv}</td>;
  if (hasFailed || !values || values.length === 0) {
    return <LeftAlignedDetailsCell>-</LeftAlignedDetailsCell>;
  }

  return <AsyncPopover values={values} />;
};

const INITIAL_UPLOAD_PROGRESS = -1;

const updateUploadProgress = (
  uploadProgress: DatasetUploadStatus["upload_progress"],
  lastUploadProgress: DatasetUploadStatus["upload_progress"],
  setLastUploadProgress: React.Dispatch<React.SetStateAction<number>>
) => {
  if (lastUploadProgress !== uploadProgress) {
    if (
      uploadProgress === 1.0 &&
      lastUploadProgress !== INITIAL_UPLOAD_PROGRESS
    ) {
      DatasetUploadToast.show({
        icon: IconNames.TICK,
        intent: Intent.SUCCESS,
        message:
          "Upload was successful. Your file is being processed which will continue in the background, even if you close this window.",
      });
    }
    setLastUploadProgress(uploadProgress);
  }
};

const handleFail = (
  datasetStatus: DatasetUploadStatus,
  fileName: string | undefined,
  setHasFailed: React.Dispatch<React.SetStateAction<boolean>>,
  queryCache: QueryCache
) => {
  queryCache.cancelQueries([USE_DATASET_STATUS, datasetStatus.dataset_id]);
  // If there is no filename present, we know there's been a refresh
  if (fileName)
    DatasetUploadToast.show({
      action:
        datasetStatus.validation_status === VALIDATION_STATUS.INVALID
          ? {
              href:
                "https://github.com/chanzuckerberg/cellxgene/blob/main/dev_docs/schema_guide.md",
              target: "_blank",
              text: "Learn More",
            }
          : {},
      intent: Intent.DANGER,
      message:
        datasetStatus.upload_status === UPLOAD_STATUS.FAILED
          ? "There was a problem uploading your file. Please try again."
          : "You must validate your dataset locally before uploading. We provide a local CLI script to do this.",
    });
  setHasFailed(true);
};

const renderUploadStatus = (datasetStatus: DatasetUploadStatus) => {
  if (datasetStatus.upload_status === UPLOAD_STATUS.FAILED)
    return (
      <DatasetStatusTag intent={Intent.DANGER}>
        <Icon iconSize={16} icon={IconNames.ISSUE} />
        Upload Error
      </DatasetStatusTag>
    );
  if (datasetStatus.validation_status === VALIDATION_STATUS.INVALID)
    return (
      <DatasetStatusTag intent={Intent.DANGER}>
        <Icon iconSize={16} icon={IconNames.ISSUE} />
        Validation Error
      </DatasetStatusTag>
    );
  if (datasetStatus.upload_progress === 1.0)
    return (
      <DatasetStatusTag>
        <Spinner intent={Intent.PRIMARY} size={16} />
        {datasetStatus.validation_status === VALIDATION_STATUS.VALIDATING &&
          "Validating..."}
      </DatasetStatusTag>
    );

  return (
    <DatasetStatusTag>
      <Spinner
        intent={Intent.PRIMARY}
        value={datasetStatus.upload_progress}
        size={16}
      />

      {`Uploading (${Math.round(datasetStatus.upload_progress * 100)}%)`}
    </DatasetStatusTag>
  );
};

const DatasetRow: FC<Props> = ({ dataset, checkHandler, file }) => {
  const {
    tissue,
    assay,
    disease,
    organism,
    cell_count,
  } = aggregateDatasetsMetadata([dataset]);
  const queryCache = useQueryCache();

  let datasetStatus = dataset.processing_status;
  const queryResult = useDatasetStatus(dataset.id);
  const [lastUploadProgress, setLastUploadProgress] = useState(
    INITIAL_UPLOAD_PROGRESS
  );
  const [hasFailed, setHasFailed] = useState(false);

  let isLoading = false;

  let name = dataset.name;
  if (!name) {
    name = file?.name ?? dataset.id;
  }

  // TODO: When checking for conversion, will have to stop polling when conversion is done and there is no need for anymore checks

  // If there is no name on the dataset the conversion and upload process hasn't completed
  // Assign a temp name and begin polling the status endpoint
  // This should be replaced with a signifier from the backend instead of relying on name population
  if (!hasFailed && !dataset.name) {
    isLoading = true;
    const { isError } = queryResult;
    if (isError) console.error(queryResult.data);
    if (!queryResult.data) return null;
    datasetStatus = queryResult.data;

    if (
      datasetStatus.upload_status === UPLOAD_STATUS.FAILED ||
      datasetStatus.validation_status === VALIDATION_STATUS.INVALID
    )
      handleFail(datasetStatus, file?.name, setHasFailed, queryCache);

    updateUploadProgress(
      datasetStatus.upload_progress,
      lastUploadProgress,
      setLastUploadProgress
    );
  }
  return (
    <StyledRow>
      <DetailsCell>
        <TitleContainer>
          <Checkbox onChange={() => checkHandler(dataset.id)} />
          <div>{name}</div>
        </TitleContainer>
        {(isLoading || hasFailed) && renderUploadStatus(datasetStatus)}
      </DetailsCell>
      {conditionalPopover(tissue, isLoading, hasFailed)}
      {conditionalPopover(assay, isLoading, hasFailed)}
      {conditionalPopover(disease, isLoading, hasFailed)}
      {conditionalPopover(organism, isLoading, hasFailed)}
      {isLoading ? (
        <td>skeletonDiv</td>
      ) : (
        <RightAlignedDetailsCell>
          {hasFailed ? "-" : cell_count}
        </RightAlignedDetailsCell>
      )}
      <RightAlignedDetailsCell>
        {!isLoading && !hasFailed && (
          <Button intent={Intent.PRIMARY} outlined text="Explore" />
        )}
      </RightAlignedDetailsCell>
    </StyledRow>
  );
};

export default DatasetRow;
