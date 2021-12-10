import DatasetsGrid from "src/components/Collection/components/CollectionDatasetsGrid/components/DatasetsGrid";
import styled from "styled-components";

export const CollectionDatasetsGrid = styled(DatasetsGrid)`
  grid-template-columns:
    minmax(520px, 6.9fr) minmax(200px, 2.5fr) minmax(160px, 2fr)
    repeat(2, minmax(104px, 1.3fr)) minmax(80px, 1fr) minmax(64px, auto); /* collection column min value at 520px to allow for more dropdown action button */
`;
