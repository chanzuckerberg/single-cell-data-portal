import styled from "@emotion/styled";
import DatasetsGrid from "src/components/Collection/components/CollectionDatasetsGrid/components/DatasetsGrid";

export const CollectionDatasetsGrid = styled(DatasetsGrid)`
  grid-template-columns: 12fr 5fr 4fr repeat(2, 3fr) 2fr auto;

  th,
  td {
    word-break: break-word; /* word break on columns; maintains grid fr allocation on small viewports */
  }
`;
