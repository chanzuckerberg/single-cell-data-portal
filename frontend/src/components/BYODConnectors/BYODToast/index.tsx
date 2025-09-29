import React, { useState, useEffect } from "react";
import { Notification } from "@czi-sds/components";
import { StyledToast, StyledMessage } from "./style";
import { StyledButton } from "../style";
import { TOAST_BUTTON_TEXT } from "../constants";
import { dismissToast, isToastDismissed } from "./utils";
import { useBYODModal } from "src/contexts/BYODModalContext";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { FEATURES } from "src/common/featureFlags/features";
import SparkleIcon from "src/common/images/sparkle-icon.svg";

const BYOD_TOAST_ID = "byod-ai-workspace";
const TOAST_MESSAGE =
  "Visualize and annotate your own single cell data on CZI Platform's AI Workspace.";

export default function BYODToast(): JSX.Element | null {
  const [isVisible, setIsVisible] = useState(false);
  const { openModal } = useBYODModal();
  const isBYODEnabled = useFeatureFlag(FEATURES.BYOD);

  useEffect(() => {
    setIsVisible(!isToastDismissed(BYOD_TOAST_ID));
  }, []);

  const handleClose = () => {
    dismissToast(BYOD_TOAST_ID);
    setIsVisible(false);
  };

  const handleLinkClick = () => {
    openModal();
  };

  if (!isBYODEnabled || !isVisible) return null;

  return (
    <StyledToast>
      <Notification
        intent="info"
        slideDirection="left"
        onClose={handleClose}
        icon={<SparkleIcon />}
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
