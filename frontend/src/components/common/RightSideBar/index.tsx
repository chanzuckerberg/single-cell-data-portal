import { ButtonIcon } from "@czi-sds/components";
import { Children, memo, ReactElement, ReactNode } from "react";
import {
  Position,
  SideBar as SideBarWrapper,
} from "src/components/common/SideBar/style";
import { RightSideBarPositioner, StyledTitle, HeaderContainer } from "./style";

export interface RightSidebarProperties {
  handleClose: () => void;
  title: string | ReactElement;
}
interface Props {
  children: ReactNode;
  width?: number;
  testId?: string;
}

// Values are in CSS vh units
const FULL_MAX_HEIGHT_VH = 100;
const UPPER_SECTION_MAX_HEIGHT_VH = 60;
const DRAWER_MAX_HEIGHT_VH = FULL_MAX_HEIGHT_VH - UPPER_SECTION_MAX_HEIGHT_VH;
const DEFAULT_WIDTH_PX = 400;

export default memo(function RightSideBar({
  children,
  width = DEFAULT_WIDTH_PX,
  testId,
  ...rest
}: Props): JSX.Element | null {
  const content = Children.map(
    children,
    (child) => child as ReactElement<RightSidebarProperties>
  );

  if (!content?.length) return null;

  const isSplit = content.length === 2;

  const [child1, child2] = content;

  return (
    <SideBarWrapper
      sideBarWidth={width}
      position={Position.RIGHT}
      data-testid={testId}
      // (thuang): Allows styling via `styled`
      {...rest}
    >
      <RightSideBarPositioner
        isExpanded
        maxHeight={isSplit ? UPPER_SECTION_MAX_HEIGHT_VH : FULL_MAX_HEIGHT_VH}
      >
        <HeaderContainer>
          <StyledTitle data-testid="right-sidebar-title">
            {child1.props.title}
          </StyledTitle>
          {child1.props.handleClose && (
            <ButtonIcon
              sdsIcon="xMark"
              sdsSize="medium"
              onClick={() => child1.props.handleClose()}
              sdsType="tertiary"
              data-testid="right-sidebar-close-button"
            />
          )}
        </HeaderContainer>
        {child1}
      </RightSideBarPositioner>

      {isSplit && (
        <RightSideBarPositioner isExpanded maxHeight={DRAWER_MAX_HEIGHT_VH}>
          <HeaderContainer>
            <StyledTitle data-testid="gene-info-title-split">
              {child2.props.title}
            </StyledTitle>
            {child2.props.handleClose && (
              <ButtonIcon
                sdsIcon="xMark"
                sdsSize="medium"
                onClick={() => child2.props.handleClose()}
                sdsType="tertiary"
                data-testid="gene-info-close-button-split"
              />
            )}
          </HeaderContainer>
          {child2}
        </RightSideBarPositioner>
      )}
    </SideBarWrapper>
  );
});
