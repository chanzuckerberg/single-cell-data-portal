import { Button, Classes, Intent } from "@blueprintjs/core";
import Image from "next/image";
import DownloadSVG from "src/common/images/download.svg";

function DownloadButton({ ...props }): JSX.Element {
  return (
    <Button
      minimal
      intent={Intent.PRIMARY}
      icon={
        <span className={Classes.ICON}>
          <Image src={DownloadSVG} alt="Download" />
        </span>
      }
      {...props}
    ></Button>
  );
}

export default DownloadButton;
