import styled from "@emotion/styled";
import { fontBodyS } from "@czi-sds/components";
import {
  gray300,
  gray500,
  fontWeightSemibold,
  spacesM,
  spacesS,
  spacesXxxs,
} from "src/common/theme";

export const Grid = styled.table`
  display: grid;
  grid-auto-rows: auto;
  grid-gap: 0 ${spacesM}px;
  margin: 0;

  thead,
  tbody,
  tr {
    display: contents; /* required; allows grandchildren of grid template to appear as though direct child */
  }

  /* row lines; span across grid gap */

  tr::after {
    box-shadow: inset 0 -0.5px 0 ${gray300};
    content: "";
    height: 0.5px;
    grid-column: 1 / -1; /* spans grid column's entire set out */
    margin-top: -0.5px; /* positions box shadow 0.5px above lower bounds of tr */
  }

  /* basic head and cell styles */

  th,
  td {
    ${fontBodyS}
    border: none;
    font-feature-settings: normal; /* required; overrides layout.css specification */
    letter-spacing: -0.006em;
  }

  th {
    align-self: center;
    color: ${gray500};
    font-weight: ${fontWeightSemibold};
    margin-bottom: ${spacesS}px;
    padding: ${spacesXxxs}px 0;
  }

  td {
    padding: ${spacesM}px 0;
  }
`;
