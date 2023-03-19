import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { ReactNode } from "react";
import {
  Position,
  SideBar as SideBarWrapper,
  SideBarOpenButtonWrapper,
  SideBarPositioner,
} from "src/components/common/SideBar/style";

export const EXPANDED_WIDTH_PX = 240;

/**
 * Function prop called on toggle of side bar expanded state, if specified.
 */
export type SideBarToggleFn = (expanded: boolean) => void;

export interface Props {
  children: ReactNode;
  label: ReactNode;
  width?: number;
  position?: typeof Position[keyof typeof Position];
  SideBarWrapperComponent?: typeof SideBarWrapper;
  SideBarPositionerComponent?: typeof SideBarPositioner;
  SideBarOpenButtonWrapperComponent?: typeof SideBarOpenButtonWrapper;
  testId?: string;
  disabled?: boolean;
  handleClose: () => void;
}

export default function SideBar({
  children: content,
  label,
  width = EXPANDED_WIDTH_PX,
  position = Position.LEFT,
  SideBarWrapperComponent = SideBarWrapper,
  SideBarPositionerComponent = SideBarPositioner,
  SideBarOpenButtonWrapperComponent = SideBarOpenButtonWrapper,
  testId,
  disabled,
  handleClose,
}: Props): JSX.Element {
  const sideBarWidth = width;
  const SideBarToggleButtonWrapper = SideBarOpenButtonWrapperComponent;
  const rightIcon = IconNames.CROSS;

  return (
    <SideBarWrapperComponent
      sideBarWidth={sideBarWidth}
      position={position}
      data-test-id={testId}
    >
      <SideBarPositionerComponent isExpanded={true}>
        <SideBarToggleButtonWrapper>
          <Button
            minimal
            onClick={() => handleClose()}
            rightIcon={rightIcon}
            text={label}
            disabled={disabled}
          />
        </SideBarToggleButtonWrapper>
        {content}
      </SideBarPositionerComponent>
    </SideBarWrapperComponent>
  );
}
