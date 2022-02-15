import { Callout, Classes } from "@blueprintjs/core";
import {
  BLUE,
  GRAY,
  PT_GRID_SIZE_PX,
  PT_TEXT_COLOR,
} from "src/components/common/theme";
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

// microsoft/typescript-styled-plugin #110
const hoverSelector = `.${Classes.TAB}:not([aria-disabled="true"]):hover`;

export const TabWrapper = styled.div`
  grid-column: 1 / span 8;

  & .${Classes.TABS} {
    min-width: fit-content;
    width: 100%;
  }
  & .${Classes.TAB_LIST} {
    box-shadow: inset 0px -1px 0px rgba(16, 22, 26, 0.15);
    padding-bottom: 6px;
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

  & .${Classes.TAB_INDICATOR} {
    bottom: -6px;
  }
`;

export const StyledCallout = styled(Callout)`
  color: ${BLUE.A};
  grid-column: 1/4;
  width: fit-content;
  padding: ${PT_GRID_SIZE_PX}px ${PT_GRID_SIZE_PX * 1.5}px;
  margin-bottom: ${PT_GRID_SIZE_PX * 2}px;
`;

export const ViewCollection = styled.div`
  padding: 40px;
`;

export const CollectionHero = styled.div`
  align-items: flex-start; /* top aligns collection name with action buttons */
  column-gap: 40px;
  display: flex;

  /* collection name */
  h3 {
    color: ${PT_TEXT_COLOR};
    flex: 1;
    letter-spacing: -0.23px;
    margin: 2.5px 0 0; /* facilitates the center alignment of single-line collection name and top alignment of the first line of a multi-line collection name (line height at 25px) with action buttons (height 30px) */
  }
`;

export const CollectionDetail = styled.div`
  align-items: flex-start;
  display: grid;
  grid-template-areas: "description . metadata"; /* grid areas for collection description and metadata with a null cell token (unnamed area) for the allocation of a gutter between the two columns */
  grid-template-columns: 8fr 1fr 8fr; /* grid columns for collection description and metadata (with 1fr allocated to column separation) */
  margin: 16px 0 44px;
`;

export const CollectionDescription = styled.div`
  color: ${PT_TEXT_COLOR};
  grid-area: description;
  letter-spacing: -0.1px;
  line-height: 18px;
`;
