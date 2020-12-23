import { Classes } from "@blueprintjs/core";
import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const UploadStatusContainer = styled.div`
  color: ${GRAY.A} !important;
  display: flex;
  flex-direction: row;
  & > ${Classes.SPINNER} {
    margin-right: ${PT_GRID_SIZE_PX}px;
  }
`;
