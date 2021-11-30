import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { FC, useState } from "react";
import {
  SideBar as SideBarWrapper,
  SideBarClosedButtonWrapper,
  SideBarOpenButtonWrapper,
  SideBarPositioner,
} from "src/components/common/SideBar/style";

export interface Props {
  label: string;
  isOpen?: boolean;
  width?: number;
}

const SideBar: FC<Props> = ({
  children: content,
  label,
  isOpen = false,
  width = 240, // default side bar width 240px
}) => {
  const [isExpanded, setIsExpanded] = useState<boolean>(isOpen);
  const sideBarWidth = isExpanded ? width : 36;
  const SideBarToggleButtonWrapper = isExpanded
    ? SideBarOpenButtonWrapper
    : SideBarClosedButtonWrapper;
  const rightIcon = isExpanded
    ? IconNames.CHEVRON_LEFT
    : IconNames.CHEVRON_RIGHT;
  return (
    <SideBarWrapper sideBarWidth={sideBarWidth}>
      <SideBarPositioner isExpanded={isExpanded}>
        <SideBarToggleButtonWrapper>
          <Button
            minimal
            onClick={() => setIsExpanded(!isExpanded)}
            rightIcon={rightIcon}
            text={label}
          />
        </SideBarToggleButtonWrapper>
        {isExpanded ? content : null}
      </SideBarPositioner>
    </SideBarWrapper>
  );
};

export default SideBar;
