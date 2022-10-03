import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { ReactNode, useEffect, useState } from "react";
import {
  Position,
  SideBar as SideBarWrapper,
  SideBarClosedButtonWrapper,
  SideBarOpenButtonWrapper,
  SideBarPositioner,
} from "src/components/common/SideBar/style";

const COLLAPSED_WIDTH_PX = 36;
export const EXPANDED_WIDTH_PX = 240;

/**
 * Function prop called on toggle of side bar expanded state, if specified.
 */
export type SideBarToggleFn = (expanded: boolean) => void;

export interface Props {
  children: ReactNode;
  label: ReactNode;
  isOpen?: boolean;
  onToggle?: SideBarToggleFn;
  width?: number;
  position?: typeof Position[keyof typeof Position];
  SideBarWrapperComponent?: typeof SideBarWrapper;
  SideBarPositionerComponent?: typeof SideBarPositioner;
  testId?: string;
  disabled?: boolean;
  forceToggle?: boolean;
  wmgSideBar?: boolean;
}

export default function SideBar({
  children: content,
  label,
  isOpen = false,
  onToggle,
  width = EXPANDED_WIDTH_PX,
  position = Position.LEFT,
  SideBarWrapperComponent = SideBarWrapper,
  SideBarPositionerComponent = SideBarPositioner,
  testId,
  disabled,
  forceToggle,
  wmgSideBar,
}: Props): JSX.Element {
  const [isExpanded, setIsExpanded] = useState(isOpen);
  const sideBarWidth = isExpanded ? width : COLLAPSED_WIDTH_PX;
  const SideBarToggleButtonWrapper = isExpanded
    ? SideBarOpenButtonWrapper
    : SideBarClosedButtonWrapper;
  const rightIcon = (position === Position.LEFT ? isExpanded : !isExpanded)
    ? IconNames.CHEVRON_LEFT
    : IconNames.CHEVRON_RIGHT;

  /**
   * Handle click on open/close icon; update state.
   * @param nextExpanded - Toggled expanded state of side bar.
   */
  const handleExpandedClick = (nextExpanded: boolean) => {
    setIsExpanded(nextExpanded);
    if (onToggle) {
      onToggle(nextExpanded);
    }
  };

  useEffect(() => {
    if (!(disabled && !isExpanded) && wmgSideBar)
      handleExpandedClick(!isExpanded);
  }, [forceToggle]);

  return (
    <SideBarWrapperComponent
      sideBarWidth={sideBarWidth}
      position={position}
      data-test-id={testId}
    >
      <SideBarPositionerComponent isExpanded={isExpanded}>
        <SideBarToggleButtonWrapper>
          <Button
            minimal
            onClick={() => handleExpandedClick(!isExpanded)}
            rightIcon={rightIcon}
            text={label}
            disabled={disabled}
          />
        </SideBarToggleButtonWrapper>
        {isExpanded ? content : null}
      </SideBarPositionerComponent>
    </SideBarWrapperComponent>
  );
}
