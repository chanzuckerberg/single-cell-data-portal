import React, { useState, useEffect } from "react";
import { Notification, Icon } from "@czi-sds/components";
import { StyledToast, StyledMessage } from "./style";
import { isToastDismissed, dismissToast } from "./utils";
import { StyledButton } from "../style";
import { TOAST_BUTTON_TEXT } from "../constants";

const BYOD_TOAST_ID = "byod-ai-workspace";
const TOAST_MESSAGE =
  "Visualize and annotate your own single cell data on CZI's Platform with our new AI Workspace.";

export default function BYODToast(): JSX.Element | null {
  const [isVisible, setIsVisible] = useState(false);

  useEffect(() => {
    setIsVisible(!isToastDismissed(BYOD_TOAST_ID));
  }, []);

  const handleClose = () => {
    dismissToast(BYOD_TOAST_ID);
    setIsVisible(false);
  };

  const handleLinkClick = () => {
    console.log("Open BYOD info dialog");
  };

  if (!isVisible) return null;

  return (
    <StyledToast>
      <Notification
        intent="info"
        slideDirection="left"
        onClose={handleClose}
        icon={<Icon sdsIcon="InfoCircle" sdsSize="s" sdsType="static" />}
      >
        <StyledMessage>{TOAST_MESSAGE}</StyledMessage>
        <StyledButton
          sdsType="primary"
          sdsStyle="minimal"
          onClick={handleLinkClick}
        >
          {TOAST_BUTTON_TEXT}
        </StyledButton>
      </Notification>
    </StyledToast>
  );
}
