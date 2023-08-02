import styled from "@emotion/styled";
import { Divider } from "@mui/material";
import { CommonThemeProps } from "@czi-sds/components";
import { scrollbar } from "src/components/common/Filter/components/FilterContent/components/common/style";
import { spacesS } from "src/common/theme";

export const VIEW_LIST_ITEM_HEIGHT = 32;
export const VIEW_LIST_SUBHEADER_HEIGHT = 23;

interface PanelProps {
  panelWidth: number;
}

interface ScrollProps extends CommonThemeProps {
  maxHeight: number;
  scrollable?: boolean;
}

export const CategoryViewPanel = styled.div`
  width: 360px;
`;

export const ViewPanel = styled.div<PanelProps>`
  min-width: ${(props) =>
    `${props.panelWidth}px`}; /* required; makes allowances for list item selected state font weight changes by maintaining panel min width */
`;

export const ViewPanelScroll = styled.div<ScrollProps>`
  max-height: ${({ maxHeight }) => `${maxHeight}px`};
  overflow-y: auto;
  padding-right: ${(props) =>
    props.scrollable ? `${spacesS(props)}px` : undefined};
  ${scrollbar}
`;

export const ViewDivider = styled(Divider)`
  background-color: rgba(16, 22, 26, 0.15);
  margin-left: ${spacesS}px;
  margin-right: ${spacesS}px;
`;
