import styled from "@emotion/styled";
import { CommonThemeProps } from "@czi-sds/components";
import { scrollbar } from "src/components/common/Filter/components/FilterContent/components/common/style";
import { spacesS, gray200 } from "src/common/theme";

export const VIEW_LIST_ITEM_HEIGHT = 32;
export const VIEW_LIST_SUBHEADER_HEIGHT = 23;

interface PanelProps extends CommonThemeProps {
  panelWidth: number;
}

interface ScrollProps extends CommonThemeProps {
  maxHeight: number;
  scrollable?: boolean;
}

export const CategoryViewPanel = styled.div`
  width: 360px;
  border-right: 1px solid ${gray200};
  padding-right: ${spacesS}px;
  &:last-of-type {
    border-right: none;
    padding-right: 0;
  }
`;

export const ViewPanel = styled.div<PanelProps>`
  min-width: ${(props) =>
    `${props.panelWidth}px`}; /* required; makes allowances for list item selected state font weight changes by maintaining panel min width */
  border-right: 1px solid ${(props) => gray200(props)};
  padding-right: ${(props) => spacesS(props)};
  &:last-of-type {
    border-right: none;
    padding-right: 0;
  }
`;

export const ViewPanelScroll = styled.div<ScrollProps>`
  max-height: ${({ maxHeight }) => `${maxHeight}px`};
  min-width: 250px;
  overflow-y: auto;
  padding-right: ${(props) =>
    props.scrollable ? `${spacesS(props)}px` : undefined};
  ${scrollbar}
`;
