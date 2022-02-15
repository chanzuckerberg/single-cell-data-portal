import Grid from "src/components/common/Grid";
import styled from "styled-components";

export const DatasetsGrid = styled(Grid)`
  grid-template-columns:
    minmax(0, 8fr) minmax(0, 5fr) minmax(0, 4fr)
    repeat(2, minmax(0, 3fr))
    minmax(0, 2fr) 56px;
`;
