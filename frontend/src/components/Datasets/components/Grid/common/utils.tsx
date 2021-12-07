import { IButtonProps, Intent } from "@blueprintjs/core";
import React from "react";
import ActionButton from "src/components/common/Grid/components/ActionButton";
import DatasetDownloadSvg from "src/components/Datasets/components/Grid/components/DatasetDownloadSvg";

/**
 * Returns dataset download action button.
 * @param downloadProps
 * @returns dataset download action button
 */
export function DownloadButton(downloadProps: IButtonProps): JSX.Element {
  return (
    <ActionButton
      iconSvg={<DatasetDownloadSvg />}
      intent={Intent.PRIMARY} // required for BP primary hover state
      {...downloadProps}
    />
  );
}
