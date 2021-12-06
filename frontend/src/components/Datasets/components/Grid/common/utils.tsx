import { IButtonProps } from "@blueprintjs/core";
import React from "react";
import ActionButton from "src/components/common/Grid/components/ActionButton";
import downloadSVG from "/src/common/images/download-blue.svg";

/**
 * Returns dataset download action button.
 * @param downloadProps
 * @returns dataset download action button
 */
export function DownloadButton(downloadProps: IButtonProps): JSX.Element {
  return <ActionButton imageProps={downloadSVG} {...downloadProps} />;
}
