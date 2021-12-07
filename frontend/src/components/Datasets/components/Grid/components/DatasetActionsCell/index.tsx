import { Intent, Position, Tooltip } from "@blueprintjs/core";
import React from "react";
import DownloadDataset from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DownloadDataset";
import { OVER_MAX_CELL_COUNT_TOOLTIP } from "src/components/common/Grid/common/constants";
import ActionButton from "src/components/common/Grid/components/ActionButton";
import ActionsCell from "src/components/common/Grid/components/ActionsCell";
import { DownloadButton } from "src/components/Datasets/components/Grid/common/utils";
import exploreSVG from "/src/common/images/explore-blue.svg";

interface Props {
  datasetId: string;
  isOverMaxCellCount: boolean;
  isRDSSkipped: boolean;
  name: string;
  tombstone: boolean;
  explorerUrl: string;
}

export default function DatasetsActionsCell({
  datasetId,
  explorerUrl,
  isOverMaxCellCount,
  isRDSSkipped,
  name,
  tombstone,
}: Props): JSX.Element {
  const exploreTooltipContent = isOverMaxCellCount
    ? OVER_MAX_CELL_COUNT_TOOLTIP
    : "Explore";
  const exploreTooltipIntent = isOverMaxCellCount ? Intent.DANGER : undefined;

  // Determine if Explorer link and tooltip is disabled.
  const explorerDisabled = tombstone || !explorerUrl;

  return (
    <ActionsCell>
      <DownloadDataset
        Button={DownloadButton}
        datasetId={datasetId}
        explorerUrl={explorerUrl}
        isDisabled={tombstone || !explorerUrl}
        isRDSSkipped={isRDSSkipped}
        name={name}
      />
      <Tooltip
        boundary="viewport"
        content={exploreTooltipContent}
        disabled={explorerDisabled}
        intent={exploreTooltipIntent}
        position={Position.TOP}
      >
        <ActionButton
          data-test-id="view-dataset-link"
          imageProps={exploreSVG}
          isDisabled={explorerDisabled || isOverMaxCellCount}
          // @ts-expect-error -- revisit rel typing
          rel="noopener"
          target="_blank"
          url={explorerUrl}
        />
      </Tooltip>
    </ActionsCell>
  );
}
