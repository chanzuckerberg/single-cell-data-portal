import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { ReactNode, useState } from "react";
import {
  Position,
  SideBar as SideBarWrapper,
  SideBarClosedButtonWrapper,
  SideBarOpenButtonWrapper,
  SideBarPositioner,
} from "src/components/common/SideBar/style";

const COLLAPSED_WIDTH_PX = 36;
const EXPANDED_WIDTH_PX = 240;

export interface Props {
  children: ReactNode;
  label: string;
  isOpen?: boolean;
  width?: number;
  position?: typeof Position[keyof typeof Position];
  SideBarWrapperComponent?: typeof SideBarWrapper;
  SideBarPositionerComponent?: typeof SideBarPositioner;
}

export default function SideBar({
  children: content,
  label,
  isOpen = false,
  width = EXPANDED_WIDTH_PX,
  position = Position.LEFT,
  SideBarWrapperComponent = SideBarWrapper,
  SideBarPositionerComponent = SideBarPositioner,
}: Props): JSX.Element {
  const [isExpanded, setIsExpanded] = useState(isOpen);
  const sideBarWidth = isExpanded ? width : COLLAPSED_WIDTH_PX;
  const SideBarToggleButtonWrapper = isExpanded
    ? SideBarOpenButtonWrapper
    : SideBarClosedButtonWrapper;
  const rightIcon = (position === Position.LEFT ? isExpanded : !isExpanded)
    ? IconNames.CHEVRON_LEFT
    : IconNames.CHEVRON_RIGHT;

  return (
    <SideBarWrapperComponent sideBarWidth={sideBarWidth} position={position}>
      <SideBarPositionerComponent isExpanded={isExpanded}>
        <SideBarToggleButtonWrapper>
          <Button
            minimal
            onClick={() => setIsExpanded(!isExpanded)}
            rightIcon={rightIcon}
            text={label}
          />
        </SideBarToggleButtonWrapper>
        {isExpanded ? content : null}
      </SideBarPositionerComponent>
    </SideBarWrapperComponent>
  );
}
