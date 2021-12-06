import { Intent, Position, Tooltip } from "@blueprintjs/core";
import React from "react";
import { Dataset, DatasetDeployment } from "src/common/entities";
import DownloadDataset from "src/components/Collections/components/Dataset/components/DownloadDataset";
import { OVER_MAX_CELL_COUNT_TOOLTIP } from "src/components/common/Grid/common/constants";
import {
  checkIsOverMaxCellCount,
  hasCXGFile,
} from "src/components/common/Grid/common/utils";
import ActionButton from "src/components/common/Grid/components/ActionButton";
import ActionsCell from "src/components/common/Grid/components/ActionsCell";
import { DownloadButton } from "src/components/Datasets/components/Grid/common/utils";
import exploreSVG from "/src/common/images/explore-blue.svg";

interface Props {
  cellCount: number;
  dataAssets: Dataset["dataset_assets"];
  datasetDeployments: DatasetDeployment[];
  isRDSSkipped: boolean;
  name: string;
  tombstone: boolean;
  url: string;
}

export default function DatasetsActionsCell({
  cellCount,
  dataAssets,
  datasetDeployments,
  isRDSSkipped,
  name,
  tombstone,
  url,
}: Props): JSX.Element {
  const isOverMaxCellCount = checkIsOverMaxCellCount(cellCount);
  const exploreTooltipContent = isOverMaxCellCount
    ? OVER_MAX_CELL_COUNT_TOOLTIP
    : "Explore";
  const exploreTooltipIntent = isOverMaxCellCount ? Intent.DANGER : undefined;

  return (
    <ActionsCell>
      <DownloadDataset
        Button={DownloadButton}
        dataAssets={dataAssets}
        isDisabled={tombstone}
        isRDSSkipped={isRDSSkipped}
        name={name}
      />
      {hasCXGFile(datasetDeployments) && (
        <Tooltip
          boundary="viewport"
          content={exploreTooltipContent}
          disabled={tombstone}
          intent={exploreTooltipIntent}
          position={Position.TOP}
        >
          <ActionButton
            data-test-id="view-dataset-link"
            imageProps={exploreSVG}
            isDisabled={tombstone || isOverMaxCellCount}
            // @ts-expect-error -- revisit rel typing
            rel="noopener"
            target="_blank"
            url={url}
          />
        </Tooltip>
      )}
    </ActionsCell>
  );
}
