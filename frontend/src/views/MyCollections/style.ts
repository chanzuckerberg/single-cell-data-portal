import styled from "styled-components";
import CreateCollection from "../../components/CreateCollectionModal";

export const TitleAndDescription = styled.div`
  grid-column: 1/6;
`;

export const StyledCreateCollection = styled(CreateCollection)`
  grid-column: 8;
  align-self: end;
`;
