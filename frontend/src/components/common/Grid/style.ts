import styled from "@emotion/styled";
import { CommonThemeProps, fontBodyS, getColors, getFontWeights } from "czifui";

const gray300 = (props: CommonThemeProps) => getColors(props)?.gray[300];
const gray500 = (props: CommonThemeProps) => getColors(props)?.gray[500];
const semiBold = (props: CommonThemeProps) => getFontWeights(props)?.semibold;

export const Grid = styled.table`
  display: grid;
  grid-auto-rows: auto;
  grid-gap: 0 12px;
  margin: 0;

  thead,
  tbody,
  tr {
    display: contents; /* required; allows grandchildren of grid template to appear as though direct child */
  }

  /* row lines; span across grid gap */

  tr::after {
    box-shadow: inset 0px -0.5px 0px ${gray300};
    content: "";
    height: 0.5px;
    grid-column: 1 / -1; /* spans grid column's entire set out */
    margin-top: -0.5px; /* positions box shadow 0.5px above lower bounds of tr */
  }

  /* basic head and cell styles */

  th,
  td {
    ${fontBodyS};
    border: none;
    font-feature-settings: normal; /* required; overrides layout.css specification */
    letter-spacing: -0.006em;
  }

  th {
    align-self: center;
    color: ${gray500};
    font-weight: ${semiBold};
    margin-bottom: 8px;
    padding: 2px 0;
  }

  td {
    padding: 12px 0;
  }
`;
