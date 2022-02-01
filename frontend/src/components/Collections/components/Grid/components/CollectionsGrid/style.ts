import Grid from "src/components/common/Grid";
import styled from "styled-components";

export const CollectionsGrid = styled(Grid)`
  grid-template-columns:
    minmax(0, 8fr)
    minmax(0, 5.4fr) repeat(2, minmax(0, 5fr)) minmax(0, 3fr);

  th:first-of-type,
  td:first-of-type {
    grid-column: 1 / 3; /* review grid column allocation for collection name when publications column added */
  }
`;
