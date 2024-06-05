import { ReactNode, useEffect, useState } from "react";
import { Icon } from "@czi-sds/components";
import {
  Position,
  SideBar as SideBarWrapper,
  SideBarPositioner,
  ToggleButton,
  ToggleButtonText,
} from "src/components/common/SideBar/style";
import { noop } from "src/common/constants/utils";

const COLLAPSED_WIDTH_PX = 36;
export const FILTERS_PANEL_EXPANDED_WIDTH_PX = 240;

/**
 * Function prop called on toggle of side bar expanded state, if specified.
 */
export type SideBarToggleFn = (expanded: boolean) => void;

export interface SidebarProps {
  children: ReactNode;
  className?: string;
  label: ReactNode;
  isOpen?: boolean;
  onToggle?: SideBarToggleFn;
  width?: number;
  position?: (typeof Position)[keyof typeof Position];
  SideBarWrapperComponent?: typeof SideBarWrapper;
  SideBarPositionerComponent?: typeof SideBarPositioner;
  testId?: string;
  disabled?: boolean;
  wmgSideBar?: boolean;
  truncatedLabel?: string;
  onWidthChange?: (width: number) => void;
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
  wmgSideBar,
  truncatedLabel = "",
  onWidthChange = noop,
}: SidebarProps): JSX.Element {
  // seve: wmgSideBar does not have isOpen prop, so we need to set default to true/open
  const [isExpanded, setIsExpanded] = useState(isOpen || !!wmgSideBar);
  const sideBarWidth = isExpanded ? width : COLLAPSED_WIDTH_PX;

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
    onWidthChange(sideBarWidth);
  }, [sideBarWidth, onWidthChange]);

  return (
    <SideBarWrapperComponent
      className={className}
      sideBarWidth={sideBarWidth}
      position={position}
      data-testid={testId}
    >
      <SideBarPositionerComponent isExpanded={isExpanded}>
        <ToggleButton
          data-testid="side-bar-toggle-button"
          disabled={disabled}
          endIcon={
            <Icon
              color="gray"
              sdsIcon={
                (position === Position.LEFT ? isExpanded : !isExpanded)
                  ? "ChevronLeft"
                  : "ChevronRight"
              }
              sdsSize="s"
              sdsType="static"
              shade={500}
            />
          }
          fullWidth
          isAllCaps={false}
          isExpanded={isExpanded}
          onClick={() => handleExpandedClick(!isExpanded)}
          sdsStyle="minimal"
          sdsType="secondary"
        >
          <ToggleButtonText>
            {!isExpanded && truncatedLabel ? truncatedLabel : label}
          </ToggleButtonText>
        </ToggleButton>
        {isExpanded ? content : null}
      </SideBarPositionerComponent>
    </SideBarWrapperComponent>
  );
}
