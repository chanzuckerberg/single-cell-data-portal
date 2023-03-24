import { ButtonIcon } from "czifui";
import { Children, memo, ReactElement, ReactNode } from "react";
import {
  Position,
  SideBar as SideBarWrapper,
} from "src/components/common/SideBar/style";
import { CELL_INFO_SIDEBAR_WIDTH_PX } from "../CellInfoSideBar/style";
import { RightSideBarPositioner, StyledTitle, HeaderContainer } from "./style";

export interface RightSidebarProperties {
  handleClose: () => void;
  title: string;
}
interface Props {
  children: ReactNode;
  width?: number;
  testId?: string;
}

// Values are in CSS vh units
const FULL_MAX_HEIGHT = 100;
const UPPER_SECTION_MAX_HEIGHT = 60;
const DRAWER_MAX_HEIGHT = FULL_MAX_HEIGHT - UPPER_SECTION_MAX_HEIGHT;

export default memo(function RightSideBar({
  children,
  width = CELL_INFO_SIDEBAR_WIDTH_PX,
  testId,
}: Props): JSX.Element | null {
  const content = Children.map(
    children,
    (child) => child as ReactElement<RightSidebarProperties>
  );

  if (!content) return null;

  const isSplit = content.length === 2;

  return (
    <SideBarWrapper
      sideBarWidth={width}
      position={Position.RIGHT}
      data-test-id={testId}
    >
      <RightSideBarPositioner
        isExpanded={true}
        maxHeight={isSplit ? UPPER_SECTION_MAX_HEIGHT : FULL_MAX_HEIGHT}
      >
        <HeaderContainer>
          <StyledTitle data-test-id="right-sidebar-title">
            {content[0].props.title}
          </StyledTitle>
          <ButtonIcon
            sdsIcon="xMark"
            sdsSize="medium"
            onClick={() => content[0].props.handleClose()}
            sdsType="tertiary"
            data-test-id="right-sidebar-close-button"
          />
        </HeaderContainer>
        {content[0]}
      </RightSideBarPositioner>

      {isSplit && (
        <RightSideBarPositioner isExpanded={true} maxHeight={DRAWER_MAX_HEIGHT}>
          <HeaderContainer>
            <StyledTitle data-test-id="gene-info-title-split">
              {content[1].props.title}
            </StyledTitle>
            <ButtonIcon
              sdsIcon="xMark"
              sdsSize="medium"
              onClick={() => content[1].props.handleClose()}
              sdsType="tertiary"
              data-test-id="gene-info-close-button-split"
            />
          </HeaderContainer>
          {content[1]}
        </RightSideBarPositioner>
      )}
    </SideBarWrapper>
  );
});
