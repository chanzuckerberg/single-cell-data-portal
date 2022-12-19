import styled from "@emotion/styled";
import { Divider } from "@mui/material";
import { CommonThemeProps, ListSubheader } from "czifui";
import { scrollbar } from "src/components/common/Filter/common/style";
import { GRAY } from "src/components/common/theme";

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
  padding-right: ${({ scrollable }) => (scrollable ? "8px" : undefined)};
  ${scrollbar}
`;

export const ViewDivider = styled(Divider)`
  background-color: rgba(16, 22, 26, 0.15);
  margin-left: 8px;
  margin-right: 8px;
`;

export const ViewHeader = styled(ListSubheader)`
  /* TODO(cc) remove && after upgrading SDS version to have this commit https://github.com/chanzuckerberg/sci-components/pull/201 */
  && {
    background-color: #ffffff;
    color: ${GRAY.A};
    cursor: default;
    font-size: 11px;
    font-weight: 500;
    letter-spacing: 0.03em;
    line-height: 15px;
    margin-bottom: 0;
    padding: 0 0 8px 8px;
    text-transform: uppercase;
  }
`;
