import { Alert, Classes } from "@blueprintjs/core";
import { BLUE, PT_GRID_SIZE_PX, RED } from "src/components/common/theme";
import styled from "styled-components";

export const CollectionInfo = styled.div`
  grid-column: 1 / span 5;
  margin-bottom: ${3 * PT_GRID_SIZE_PX}px;
`;

export const Description = styled.div`
  font-size: 14px;
  font-style: normal;
  font-weight: 400;
  line-height: 18px;
  letter-spacing: -0.1px;
  text-align: left;
  width: 90ch;
`;

export const LinkContainer = styled.div`
  display: grid;
  grid-template-columns: max-content auto;
  column-gap: ${2 * PT_GRID_SIZE_PX}px;
  row-gap: ${PT_GRID_SIZE_PX}px;
  margin-top: ${PT_GRID_SIZE_PX}px;
`;

export const DatasetContainer = styled.div`
  grid-column: 1 / span 8;
`;

export const StyledDiv = styled.div`
  display: flex;
  flex-direction: row;
  & > :not(:last-child) {
    margin-right: ${PT_GRID_SIZE_PX}px;
  }
`;

export const StyledAlert = styled(Alert)`
  .${Classes.BUTTON} {
    background: ${RED.C};
    box-shadow: none !important;
    :not(.${Classes.INTENT_DANGER}) {
      color: ${BLUE.C};
      background: none;
    }
  }
`;
