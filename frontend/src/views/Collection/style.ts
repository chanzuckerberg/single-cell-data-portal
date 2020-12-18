import { PT_GRID_SIZE_PX } from "src/components/common/theme";
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
  border: 1px solid #e1e8ed;
  border-radius: 3px;
`;
