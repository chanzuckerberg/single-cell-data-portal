import { Button, Classes, IButtonProps, Intent } from "@blueprintjs/core";
import Image from "next/image";
import DownloadSVG from "src/common/images/download.svg";

function DownloadButton(props: IButtonProps): JSX.Element {
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
