import { Classes, Radio } from "@blueprintjs/core";
import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const UploadStatusContainer = styled.div`
  color: ${GRAY.A} !important;
`;

export const DatasetTitleCell = styled.td`
  vertical-align: middle !important;
`;

export const TitleContainer = styled.div`
  display: flex;
  flex-direction: row;
  & > .${Classes.CHECKBOX} {
    margin: 0 ${PT_GRID_SIZE_PX}px 0 0;
  }
  white-space: normal;
  align-items: center;
  & > .${Classes.POPOVER_WRAPPER} {
    margin-left: ${PT_GRID_SIZE_PX}px;
  }
  margin: ${2 * PT_GRID_SIZE_PX}px 0;
`;

export const StyledRadio = styled(Radio)`
  margin: 0;
`;
