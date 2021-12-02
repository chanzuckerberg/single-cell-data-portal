import Grid from "src/components/common/Grid";
import styled from "styled-components";

export const CollectionsGrid = styled(Grid)`
  grid-template-columns: minmax(552px, 4.6fr) repeat(2, minmax(200px, 1.5fr)) minmax(
      120px,
      1fr
    ); /* TODO(cc) review grid template; publications omitted for now */
`;
