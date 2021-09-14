import { Button, Classes, Intent } from "@blueprintjs/core";
import DownloadSVG from "src/common/images/download.svg";

function DownloadButton({ ...props }) {
  return (
    <Button
      minimal
      intent={Intent.PRIMARY}
      icon={
        <span className={Classes.ICON}>
          <img src={String(DownloadSVG)} alt="Download" />
        </span>
      }
      {...props}
    ></Button>
  );
}

export default DownloadButton;
