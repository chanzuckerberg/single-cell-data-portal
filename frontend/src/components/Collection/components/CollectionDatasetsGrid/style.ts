import styled from "@emotion/styled";
import DatasetsGrid from "src/components/Collection/components/CollectionDatasetsGrid/components/DatasetsGrid";
import { CommonThemeProps } from "@czi-sds/components";

const DEFAULT_GRID_TEMPLATE_COLUMNS = "12fr 5fr 4fr repeat(2, 3fr) 2fr auto";
const REORDER_GRID_TEMPLATE_COLUMNS =
  "auto 12fr 5fr 4fr repeat(2, 3fr) 2fr auto";

interface GridProps extends CommonThemeProps {
  isReorder?: boolean;
}

export const CollectionDatasetsGrid = styled(DatasetsGrid)<GridProps>`
  grid-template-columns: ${({ isReorder }) =>
    isReorder ? REORDER_GRID_TEMPLATE_COLUMNS : DEFAULT_GRID_TEMPLATE_COLUMNS};

  th,
  td {
    word-break: break-word; /* word break on columns; maintains grid fr allocation on small viewports */
  }
`;
