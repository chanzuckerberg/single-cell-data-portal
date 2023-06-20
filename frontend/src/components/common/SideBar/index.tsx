import { ReactNode, useEffect, useState } from "react";
import { Button, Icon } from "@czi-sds/components";
import {
  Position,
  SideBar as SideBarWrapper,
  SideBarClosedButtonWrapper,
  SideBarOpenButtonWrapper,
  SideBarPositioner,
  ToggleButtonText,
} from "src/components/common/SideBar/style";

const COLLAPSED_WIDTH_PX = 36;
export const FILTERS_PANEL_EXPANDED_WIDTH_PX = 240;

/**
 * Function prop called on toggle of side bar expanded state, if specified.
 */
export type SideBarToggleFn = (expanded: boolean) => void;

export interface Props {
  children: ReactNode;
  className?: string;
  label: ReactNode;
  isOpen?: boolean;
  onToggle?: SideBarToggleFn;
  width?: number;
  position?: typeof Position[keyof typeof Position];
  SideBarWrapperComponent?: typeof SideBarWrapper;
  SideBarPositionerComponent?: typeof SideBarPositioner;
  testId?: string;
  disabled?: boolean;
  forceOpen?: boolean;
  wmgSideBar?: boolean;
  truncatedLabel?: string;
}

export default function SideBar({
  children: content,
  className,
  label,
  isOpen = false,
  onToggle,
  width = FILTERS_PANEL_EXPANDED_WIDTH_PX,
  position = Position.LEFT,
  SideBarWrapperComponent = SideBarWrapper,
  SideBarPositionerComponent = SideBarPositioner,
  testId,
  disabled,
  forceOpen,
  wmgSideBar,
  truncatedLabel = "",
}: Props): JSX.Element {
  const [isExpanded, setIsExpanded] = useState(isOpen);
  const sideBarWidth = isExpanded ? width : COLLAPSED_WIDTH_PX;
  const SideBarToggleButtonWrapper = isExpanded
    ? SideBarOpenButtonWrapper
    : SideBarClosedButtonWrapper;

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
    if (!(disabled && !isExpanded) && wmgSideBar) handleExpandedClick(true);
  }, [forceOpen]);

  return (
    <SideBarWrapperComponent
      className={className}
      sideBarWidth={sideBarWidth}
      position={position}
      data-testid={testId}
    >
      <SideBarPositionerComponent isExpanded={isExpanded}>
        <SideBarToggleButtonWrapper>
          <Button
            data-testid="side-bar-toggle-button"
            disabled={disabled}
            endIcon={
              <Icon
                sdsIcon={
                  (position === Position.LEFT ? isExpanded : !isExpanded)
                    ? "chevronLeft"
                    : "chevronRight"
                }
                sdsSize="l"
                sdsType="button"
              />
            }
            onClick={() => handleExpandedClick(!isExpanded)}
            sdsStyle="minimal"
            sdsType="minimal"
            size="large"
          >
            <ToggleButtonText>
              {!isExpanded && truncatedLabel ? truncatedLabel : label}
            </ToggleButtonText>
          </Button>
        </SideBarToggleButtonWrapper>
        {isExpanded ? content : null}
      </SideBarPositionerComponent>
    </SideBarWrapperComponent>
  );
}
