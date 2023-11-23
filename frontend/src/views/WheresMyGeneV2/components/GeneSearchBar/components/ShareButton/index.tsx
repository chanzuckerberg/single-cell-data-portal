// import { Notification } from "@czi-sds/components";
// import { noop } from "src/common/constants/utils";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../SaveExport";
import {
  StyledButtonDiv,
  StyledButtonIcon,
  // StyledIcon,
  StyledLabel,
  // StyledNotificationDetails,
  // StyledNotificationLabel,
  // StyledNotificationWrapper,
} from "./style";
import { useConnect } from "./connect";

export default function ShareButton(): JSX.Element {
  const { copyShareUrl, selectedGenes } = useConnect();

  return (
    // <>
    //   {showURLCopyNotification > 0 && (
    //     <StyledNotificationWrapper>
    //       <Notification
    //         key={showURLCopyNotification}
    //         autoDismiss={5000}
    //         onClose={noop}
    //         slideDirection={"left"}
    //         intent={"info"}
    //         icon={
    //           <StyledIcon sdsIcon={"link"} sdsSize={"s"} sdsType={"static"} />
    //         }
    //       >
    //         <StyledNotificationLabel data-testid="share-link-notification">
    //           Share link copied
    //         </StyledNotificationLabel>
    //         <StyledNotificationDetails>
    //           We regularly expand our single cell data corpus to improve
    //           results. Downloaded data and figures may differ in the future.
    //         </StyledNotificationDetails>
    //       </Notification>
    //     </StyledNotificationWrapper>
    //   )}
    <StyledButtonDiv className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
      <StyledLabel>Share</StyledLabel>
      <StyledButtonIcon
        data-testid={"share-button"}
        onClick={copyShareUrl}
        sdsSize="medium"
        sdsType="primary"
        sdsIcon="share"
        disabled={selectedGenes.length === 0}
      />
    </StyledButtonDiv>
    // </>
  );
}
