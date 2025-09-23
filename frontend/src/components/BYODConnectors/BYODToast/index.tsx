import React, { useState, useEffect } from "react";
import { Notification, Icon } from "@czi-sds/components";
import { StyledToast, StyledMessage } from "./style";
import { isToastDismissed, dismissToast } from "./utils";
import { StyledButton } from "../style";

const BYOD_TOAST_ID = "byod-ai-workspace";
const BYOD_MESSAGE =
  "Visualize and annotate your own single cell data on CZI's Platform with our new AI Workspace.";
const BYOD_LINK_TEXT = "Learn more about AI Workspace";
const BYOD_LINK_URL = "https://cellxgene.cziscience.com/docs/";

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
    window.open(BYOD_LINK_URL, "_blank", "noopener,noreferrer");
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
        <StyledMessage>{BYOD_MESSAGE}</StyledMessage>
        <StyledButton
          sdsType="primary"
          sdsStyle="minimal"
          onClick={handleLinkClick}
        >
          {BYOD_LINK_TEXT}
        </StyledButton>
      </Notification>
    </StyledToast>
  );
}
