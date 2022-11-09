import { IButtonProps, Intent } from "@blueprintjs/core";
import React from "react";
import { StyledOutlineButton } from "src/components/common/Button/common/style";

/**
 * Returns dataset download action button.
 * @param downloadProps
 * @returns dataset download action button
 */
export function DownloadButton(downloadProps: IButtonProps): JSX.Element {
  return (
    <StyledOutlineButton
      intent={Intent.PRIMARY}
      outlined
      text="Download"
      {...downloadProps}
    />
  );
}
