import { Button, Intent } from "@blueprintjs/core";
import React from "react";
import DownloadSVG from "src/common/images/download.svg";

function DownloadButton({ ...props }) {
  return (
    <Button
      minimal
      intent={Intent.PRIMARY}
      icon={<img alt="Download" src={String(DownloadSVG)} />}
      {...props}
    ></Button>
  );
}

export default DownloadButton;
