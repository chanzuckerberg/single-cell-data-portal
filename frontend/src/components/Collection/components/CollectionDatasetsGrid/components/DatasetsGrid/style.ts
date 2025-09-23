import styled from "@emotion/styled";
import { Grid } from "src/components/common/Grid/style";
import { CommonThemeProps } from "@czi-sds/components";
import { css, SerializedStyles } from "@emotion/react";

interface DatasetsGridProps extends CommonThemeProps {
  dragAndDropStyles?: SerializedStyles;
}

export const DatasetsGrid = styled(Grid)<DatasetsGridProps>`
  ${({ dragAndDropStyles }) =>
    dragAndDropStyles &&
    css`
      tbody tr {
        ${dragAndDropStyles}
      }
    `}
`;
