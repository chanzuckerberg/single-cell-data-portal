import Grid from "src/components/common/Grid";
import styled from "styled-components";

export const CollectionsGrid = styled(Grid)`
  grid-template-columns:
    minmax(0, 8fr)
    minmax(0, 6fr) repeat(2, minmax(0, 5fr)) minmax(0, 3fr);
`;
