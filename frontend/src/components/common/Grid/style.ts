import styled from "@emotion/styled";
import { GRAY, PT_TEXT_COLOR } from "src/components/common/theme";

export const Grid = styled.table`
  display: grid;
  grid-auto-rows: auto;
  grid-gap: 0 16px;
  margin: 0;

  thead,
  tbody,
  tr {
    display: contents; /* required; allows grandchildren of grid template to appear as though direct child */
  }

  /* row lines; span across grid gap */
  tr::after {
    box-shadow: inset 0px -1px 0px rgba(16, 22, 26, 0.15);
    content: "";
    height: 1px;
    grid-column: 1 / -1; /* spans grid column's entire set out */
    margin-top: -1px; /* positions box shadow 1px above lower bounds of tr */
  }

  /* basic head and cell styles */
  th,
  td {
    border: none;
    -moz-font-feature-settings: normal;
    -webkit-font-feature-settings: normal;
    font-feature-settings: normal; /* required; overrides layout.css specification */
    font-size: 14px;
    letter-spacing: -0.1px;
    line-height: 18px;
  }

  th {
    color: ${GRAY.A};
    font-weight: 500;
    line-height: 20px;
    padding: 0 0 14px 0;
  }

  td {
    color: ${PT_TEXT_COLOR};
    padding: 12px 0;
  }
`;
