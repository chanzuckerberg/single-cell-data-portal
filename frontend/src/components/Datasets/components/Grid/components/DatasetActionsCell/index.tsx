import { Intent, Position, Tooltip } from "@blueprintjs/core";
import DownloadDataset from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DownloadDataset";
import { OVER_MAX_CELL_COUNT_TOOLTIP } from "src/components/common/Grid/common/constants";
import ActionsCell from "src/components/common/Grid/components/ActionsCell";
import { StyledPrimaryAnchorButton } from "src/components/common/Button/common/style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import DownloadButton from "src/components/common/Grid/components/DownloadButton";

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
        content={OVER_MAX_CELL_COUNT_TOOLTIP}
        disabled={!isOverMaxCellCount}
        intent={Intent.DANGER}
        position={Position.TOP}
      >
        <StyledPrimaryAnchorButton
          data-testid="view-dataset-link"
          disabled={explorerDisabled || isOverMaxCellCount} // TODO(cc) isOverMaxCellCount needed to be a boolean...
          href={explorerUrl}
          intent={Intent.PRIMARY}
          onClick={() =>
            track(EVENTS.DATASET_EXPLORE_CLICKED, {
              dataset_name: name,
            })
          }
          rel="noopener"
          target="_blank"
          text="Explore"
        />
      </Tooltip>
    </ActionsCell>
  );
}
