import { Classes } from "@blueprintjs/core";
import { BLUE, GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
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

export const CollectionButtons = styled.div`
  grid-column: 6 / span 3;
  justify-self: end;
  display: flex;
  flex-direction: row;
  & > :not(:last-child) {
    margin-right: ${PT_GRID_SIZE_PX * 2}px;
  }
`;

// microsoft/typescript-styled-plugin #110
const hoverSelector = `.${Classes.TAB}:not([aria-disabled="true"]):hover`;

export const TabWrapper = styled.div`
  grid-column: 1 / span 8;

  & .${Classes.TABS} {
    width: 100%;
  }
  & .${Classes.TAB_LIST} {
    box-shadow: inset 0px -1px 0px rgba(16, 22, 26, 0.15);
  }
  & .${Classes.TAB} {
    color: ${GRAY.A};
    font-style: normal;
    font-weight: normal;
    font-size: 14px;
    line-height: 18px;
    letter-spacing: -0.1px;
  }

  & .${Classes.TAB}[aria-selected="true"], ${hoverSelector} {
    color: ${BLUE.C};
  }
`;
