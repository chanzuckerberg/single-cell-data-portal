import { Intent, Position, Tooltip } from "@blueprintjs/core";
import DownloadDataset from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DownloadDataset";
import { OVER_MAX_CELL_COUNT_TOOLTIP } from "src/components/common/Grid/common/constants";
import ActionButton from "src/components/common/Grid/components/ActionButton";
import ActionsCell from "src/components/common/Grid/components/ActionsCell";
import { DownloadButton } from "src/components/Datasets/components/Grid/common/utils";
import DatasetExploreSvg from "src/components/Datasets/components/Grid/components/DatasetExploreSvg";

interface Props {
  datasetId: string;
  isOverMaxCellCount: boolean;
  name: string;
  tombstone: boolean;
  explorerUrl: string;
}

export default function DatasetsActionsCell({
  datasetId,
  explorerUrl,
  isOverMaxCellCount = false, // TODO(cc) either set isOverMaxCellCount here to false when undefined or in parent...
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
        isDisabled={tombstone || !explorerUrl}
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
          // @ts-expect-error -- revisit disabled typing
          disabled={explorerDisabled || isOverMaxCellCount} // TODO(cc) isOverMaxCellCount needed to be a boolean...
          href={explorerUrl}
          iconSvg={<DatasetExploreSvg />}
          isAnchorButton
          rel="noopener"
          target="_blank"
        />
      </Tooltip>
    </ActionsCell>
  );
}
