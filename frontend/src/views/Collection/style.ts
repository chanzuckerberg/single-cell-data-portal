import { Button } from "@blueprintjs/core";
import styled from "styled-components";

export const StyledButton = styled(Button)`
  margin-left: 16px;
`;

export const TopBar = styled.div`
  display: flex;
  width: 100%;
  justify-content: space-between;
  flex-direction: row;
`;

export const ViewGrid = styled.div`
  display: grid;
  grid-template-columns: repeat(8, 1fr);
  width: 100%;
  margin-top: 128px;
`;

export const CollectionInfo = styled.div`
  grid-column: 1 / span 5;
  margin-bottom: 24px;
`;

export const LinkContainer = styled.div`
  display: grid;
  grid-template-columns: max-content auto;
  column-gap: 16px;
  row-gap: 8px;
  margin-top: 8px;
`;

export const CollectionButtons = styled.div`
  grid-column: 6 / span 3;
  justify-self: end;
`;

export const Description = styled.div`
  font-size: 14px;
  font-style: normal;
  font-weight: 400;
  line-height: 18px;
  letter-spacing: -0.1px;
  text-align: left;
`;

export const DatasetContainer = styled.div`
  grid-column: 1 / span 8;
  border: 1px solid #e1e8ed;
  border-radius: 3px;
`;

export const CenterAlignedDiv = styled.div`
  display: flex;
  flex-direction: column;
  width: 340px;
  justify-content: center;
  align-items: center;
  margin: auto auto;
  height: 229px;
  font-weight: 400;
  line-height: 18px;
  letter-spacing: -0.10000000149011612px;
  text-align: left;
`;
