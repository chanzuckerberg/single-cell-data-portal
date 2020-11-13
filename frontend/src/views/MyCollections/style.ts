import { HTMLTable } from "@blueprintjs/core";
import styled from "styled-components";
import CreateCollection from "../../components/CreateCollectionModal";

export const TitleAndDescription = styled.div`
  grid-column: 1/6;
`;

export const StyledCreateCollection = styled(CreateCollection)`
  grid-column: 8;
  align-self: end;
`;

export const CollectionsGrid = styled(HTMLTable)`
  grid-column: 1/9;
  margin-top: 16px;
`;

export const CollectionDataCell = styled.td`
  width: calc(3 / 8 * 100%);
  text-align: left;
`;
export const CollectionHeaderCell = styled.th`
  width: calc(3 / 8 * 100%);
  text-align: left;
`;

export const LeftAlignedDataCell = styled.td`
  width: calc(1 / 8 * 100%);
  text-align: left;
`;

export const LeftAlignedHeaderCell = styled.th`
  width: calc(1 / 8 * 100%);
  text-align: left;
  margin-left: 16px;
`;

export const RightAlignedDataCell = styled.td`
  width: calc(1 / 8 * 100%);
  text-align: right;
`;

export const RightAlignedHeaderCell = styled.th`
  width: calc(1 / 8 * 100%);
  text-align: right;
  margin-left: 16px;
`;
