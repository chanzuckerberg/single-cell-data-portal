import { Classes } from "@blueprintjs/core";
import { PT_GRID_SIZE_PX, PT_TEXT_COLOR } from "src/components/common/theme";
import styled from "styled-components";
import ExplorerSvg from "./components/ExplorerSvg";

export const TitleContainer = styled.div`
  display: flex;
  flex-direction: row;
  color: ${PT_TEXT_COLOR};
  font-size: 14px;
  line-height: 18px;
  & > .${Classes.CHECKBOX} {
    margin: 0 ${PT_GRID_SIZE_PX}px 0 0;
  }
  white-space: normal;
  & > .${Classes.POPOVER_WRAPPER} {
    margin-left: ${PT_GRID_SIZE_PX}px;
  }
  padding: 0 ${PT_GRID_SIZE_PX}px;
  justify-content: space-between;
  align-items: center;
`;

export const StyledExplorerSvg = styled(ExplorerSvg)`
  margin-bottom: 4px;
`;
