import React from "react";
import { StyledOutlineButton } from "src/components/common/Button/common/style";
import { Intent } from "@blueprintjs/core";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";

interface Props {
  datasetName: string;
  onClick: () => void;
}

export default function DownloadButton({
  datasetName,
  onClick,
  ...props /* Spread props to allow for data-testid and other ButtonProps e.g. "disabled". */
}: Props): JSX.Element {
  return (
    <StyledOutlineButton
      intent={Intent.PRIMARY}
      onClick={() => {
        onClick(); // Opens download dataset modal.
        track(EVENTS.DOWNLOAD_DATA_CLICKED, {
          dataset_name: datasetName,
        }); // Tracks dataset download analytics event "DOWNLOAD_DATA_CLICKED".
      }}
      outlined
      text="Download"
      {...props}
    />
  );
}
