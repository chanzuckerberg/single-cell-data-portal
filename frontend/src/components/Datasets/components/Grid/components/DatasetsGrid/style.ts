import Grid from "src/components/common/Grid";
import styled from "styled-components";

export const DatasetsGrid = styled(Grid)`
  grid-template-columns:
    minmax(312px, 3.9fr) minmax(200px, 2.5fr) minmax(160px, 2fr)
    repeat(2, minmax(104px, 1.3fr)) minmax(80px, 1fr) 64px;
`;
