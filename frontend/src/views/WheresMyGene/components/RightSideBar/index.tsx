import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { ReactElement, ReactNode } from "react";
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
  content: { label: string; handleClose: () => void; element: ReactElement }[];
  width?: number;
  position?: typeof Position[keyof typeof Position];
  SideBarWrapperComponent?: typeof SideBarWrapper;
  SideBarPositionerComponent?: typeof SideBarPositioner;
  SideBarOpenButtonWrapperComponent?: typeof SideBarOpenButtonWrapper;
  testId?: string;
  disabled?: boolean;
}

export default function SideBar({
  content,
  width = EXPANDED_WIDTH_PX,
  position = Position.LEFT,
  SideBarWrapperComponent = SideBarWrapper,
  SideBarPositionerComponent = SideBarPositioner,
  SideBarOpenButtonWrapperComponent = SideBarOpenButtonWrapper,
  testId,
  disabled,
}: Props): JSX.Element {
  const sideBarWidth = width;
  const SideBarToggleButtonWrapper = SideBarOpenButtonWrapperComponent;
  const rightIcon = IconNames.CROSS;

  if (!content) return <></>;

  return (
    <SideBarWrapperComponent
      sideBarWidth={sideBarWidth}
      position={position}
      data-test-id={testId}
    >
      {content.length === 1 ? (
        <SideBarPositionerComponent isExpanded={true}>
          <SideBarToggleButtonWrapper>
            <Button
              minimal
              onClick={() => content[0].handleClose()}
              rightIcon={rightIcon}
              text={content[0].label}
              disabled={disabled}
            />
          </SideBarToggleButtonWrapper>
          {content[0].element}
        </SideBarPositionerComponent>
      ) : (
        <>
          <SideBarPositionerComponent
            id="c1"
            isExpanded={true}
            style={{ maxHeight: "calc(60vh - 48px)" }}
          >
            <SideBarToggleButtonWrapper>
              <Button
                minimal
                onClick={() => content[0].handleClose()}
                rightIcon={rightIcon}
                text={content[0].label}
                disabled={disabled}
              />
            </SideBarToggleButtonWrapper>
            {content[0].element}
          </SideBarPositionerComponent>

          <SideBarPositionerComponent
            id="c2"
            isExpanded={true}
            style={{
              maxHeight: "calc(40vh - 48px)",
              borderTop: "solid 1px rgb(16 22 26 / 15%)",
            }}
          >
            <SideBarToggleButtonWrapper>
              <Button
                minimal
                onClick={() => content[1].handleClose()}
                rightIcon={rightIcon}
                text={content[1].label}
                disabled={disabled}
              />
            </SideBarToggleButtonWrapper>
            {content[1].element}
          </SideBarPositionerComponent>
        </>
      )}
    </SideBarWrapperComponent>
  );
}
