import {
  Button,
  Checkbox,
  Classes,
  Colors,
  Icon,
  Intent,
  PopoverInteractionKind,
  Spinner,
  Tooltip,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import loadable from "@loadable/component";
import React, { FC, useEffect, useState } from "react";
import { QueryCache, useQueryCache } from "react-query";
import {
  CONVERSION_STATUS,
  Dataset,
  DatasetUploadStatus,
  UPLOAD_STATUS,
  VALIDATION_STATUS,
} from "src/common/entities";
import {
  useDatasetStatus,
  USE_DATASET_STATUS,
} from "src/common/queries/datasets";
import {
  DetailsCell,
  LeftAlignedDetailsCell,
  RightAlignedDetailsCell,
  StyledRow,
} from "src/components/Collections/components/Grid/common/style";
import { aggregateDatasetsMetadata } from "src/components/Collections/components/Grid/common/utils";
import { UploadingFile } from "src/components/DropboxChooser";
import DatasetUploadToast from "src/views/Collection/components/DatasetUploadToast";
import { DatasetStatusTag, StyledAnchor, TitleContainer } from "./style";

interface Props {
  dataset: Dataset;
  checkHandler: (id: string) => void;
  file?: UploadingFile;
  invalidateCollectionQuery: () => void;
}

const AsyncPopover = loadable(
  () =>
    /*webpackChunkName: 'CollectionRow/components/Popover' */ import(
      "src/components/Collections/components/Grid/components/Popover"
    )
);

const skeletonDiv = <div className={Classes.SKELETON}>PLACEHOLDER_TEXT</div>;

const conditionalPopover = (
  values: string[],
  loading: boolean,
  hasFailed: boolean
) => {
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
  datasetID: DatasetUploadStatus["dataset_id"],
  queryCache: QueryCache,
  setLastUploadProgress: React.Dispatch<React.SetStateAction<number>>,
  invalidateCollectionQuery: () => void
) => {
  if (lastUploadProgress !== uploadProgress) {
    if (uploadProgress === 1) {
      if (lastUploadProgress !== INITIAL_UPLOAD_PROGRESS) {
        DatasetUploadToast.show({
          icon: IconNames.TICK,
          intent: Intent.SUCCESS,
          message:
            "Upload was successful. Your file is being processed which will continue in the background, even if you close this window.",
        });
        invalidateCollectionQuery();
      }
      queryCache.cancelQueries([USE_DATASET_STATUS, datasetID]);
    }
    setLastUploadProgress(uploadProgress);
  }
};

type FailReturn =
  | {
      isFailed: boolean;
      error: VALIDATION_STATUS | UPLOAD_STATUS;
    }
  | {
      isFailed: false;
    };

const checkIfFailed = (datasetStatus: DatasetUploadStatus): FailReturn => {
  if (datasetStatus.validation_status === VALIDATION_STATUS.INVALID)
    return { error: VALIDATION_STATUS.INVALID, isFailed: true };
  if (datasetStatus.upload_status === UPLOAD_STATUS.FAILED)
    return { error: UPLOAD_STATUS.FAILED, isFailed: true };
  // TODO: check if conversion failed
  return { isFailed: false };
};
const checkIfLoading = (datasetStatus: DatasetUploadStatus): boolean => {
  if (checkIfFailed(datasetStatus).isFailed) return false;
  // TODO: There should be an "all done" value on datasetStatus to simplify this check
  if (
    datasetStatus.upload_status === UPLOAD_STATUS.UPLOADING ||
    datasetStatus.upload_status === UPLOAD_STATUS.WAITING
  )
    return true;
  if (datasetStatus.validation_status === VALIDATION_STATUS.VALIDATING)
    return true;
  // TODO: There should be an all encompassing conversion to simplify this part
  if (
    datasetStatus.conversion_anndata_status === CONVERSION_STATUS.CONVERTING ||
    datasetStatus.conversion_cxg_status === CONVERSION_STATUS.CONVERTING ||
    datasetStatus.conversion_rds_status === CONVERSION_STATUS.CONVERTING ||
    datasetStatus.conversion_loom_status === CONVERSION_STATUS.CONVERTING
  )
    return true;

  return false;
};

const checkIfComplete = (datasetStatus: DatasetUploadStatus): boolean => {
  return !checkIfFailed(datasetStatus) && !checkIfLoading(datasetStatus);
};

const handleFail = (
  datasetID: Dataset["id"],
  fail: FailReturn,
  setHasFailed: React.Dispatch<React.SetStateAction<FailReturn>>,
  queryCache: QueryCache
) => {
  queryCache.cancelQueries([USE_DATASET_STATUS, datasetID]);
  setHasFailed(fail);
};

const renderUploadStatus = (datasetStatus: DatasetUploadStatus) => {
  if (datasetStatus.validation_status === VALIDATION_STATUS.VALIDATING)
    return (
      <DatasetStatusTag>
        <Spinner intent={Intent.PRIMARY} size={16} />
        Validating...
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

const tooltipContent = (error: VALIDATION_STATUS | UPLOAD_STATUS) => {
  return error === VALIDATION_STATUS.INVALID ? (
    <span>
      You must validate your dataset locally before uploading. We provide a
      local CLI script to do this.{" "}
      <b>
        <StyledAnchor
          href="https://github.com/chanzuckerberg/cellxgene/blob/main/dev_docs/schema_guide.md"
          target="_blank"
        >
          Learn More
        </StyledAnchor>
      </b>
    </span>
  ) : (
    <span>There was a problem uploading your file. Please try again.</span>
  );
};

const DatasetRow: FC<Props> = ({
  dataset,
  checkHandler,
  file,
  invalidateCollectionQuery,
}) => {
  const {
    tissue,
    assay,
    disease,
    organism,
    cell_count,
  } = aggregateDatasetsMetadata([dataset]);
  const queryCache = useQueryCache();
  const queryResult = useDatasetStatus(
    dataset.id,
    checkIfLoading(dataset.processing_status)
  );
  const [lastUploadProgress, setLastUploadProgress] = useState(
    INITIAL_UPLOAD_PROGRESS
  );
  const [hasFailed, setHasFailed] = useState({
    isFailed: false,
  } as FailReturn);

  if (queryResult.isError) console.error(queryResult.error);

  const datasetStatus = queryResult.data ?? dataset.processing_status;
  const isLoading = checkIfLoading(datasetStatus);

  let name = dataset.name;
  if (name === undefined || name === "") {
    name = file?.name ?? dataset.id;
  }

  useEffect(() => {
    const failed = checkIfFailed(datasetStatus);
    if (failed.isFailed) {
      handleFail(dataset.id, failed, setHasFailed, queryCache);
      return;
    }

    // If there is no name on the dataset the conversion and upload process hasn't completed
    // Assign a temp name and begin polling the status endpoint
    // This should be replaced with a signifier from the backend instead of relying on name population
    updateUploadProgress(
      datasetStatus.upload_progress,
      lastUploadProgress,
      dataset.id,
      queryCache,
      setLastUploadProgress,
      invalidateCollectionQuery
    );
  }, [
    dataset.id,
    datasetStatus,
    invalidateCollectionQuery,
    isLoading,
    lastUploadProgress,
    queryCache,
  ]);

  return (
    <StyledRow>
      <DetailsCell>
        <TitleContainer>
          <Checkbox onChange={() => checkHandler(dataset.id)} />
          <div>{name}</div>
          {hasFailed.isFailed ? (
            <Tooltip
              intent={Intent.DANGER}
              interactionKind={PopoverInteractionKind.HOVER}
              hoverCloseDelay={500}
              content={tooltipContent(hasFailed.error)}
            >
              <Icon icon={IconNames.ISSUE} iconSize={16} color={Colors.RED3} />
            </Tooltip>
          ) : null}
        </TitleContainer>
        {isLoading && renderUploadStatus(datasetStatus)}
      </DetailsCell>
      {conditionalPopover(tissue, isLoading, hasFailed.isFailed)}
      {conditionalPopover(assay, isLoading, hasFailed.isFailed)}
      {conditionalPopover(disease, isLoading, hasFailed.isFailed)}
      {conditionalPopover(organism, isLoading, hasFailed.isFailed)}
      {isLoading ? (
        <td>{skeletonDiv}</td>
      ) : (
        <RightAlignedDetailsCell>
          {hasFailed.isFailed || !cell_count ? "-" : cell_count}
        </RightAlignedDetailsCell>
      )}
      <RightAlignedDetailsCell>
        {/* datasetStatus.conversion_cxg_status */}
        {!(isLoading || hasFailed.isFailed) && (
          <Button intent={Intent.PRIMARY} outlined text="Explore" />
        )}
      </RightAlignedDetailsCell>
    </StyledRow>
  );
};

export default DatasetRow;
