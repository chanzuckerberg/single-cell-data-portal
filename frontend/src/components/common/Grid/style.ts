import styled from "@emotion/styled";
import {
  CommonThemeProps,
  fontBodyS,
  getColors,
  getFontWeights,
  getSpaces,
} from "@czi-sds/components";

const gray300 = (props: CommonThemeProps) => getColors(props)?.gray[300];
const gray500 = (props: CommonThemeProps) => getColors(props)?.gray[500];
const semiBold = (props: CommonThemeProps) => getFontWeights(props)?.semibold;
const spacesM = (props: CommonThemeProps) => getSpaces(props)?.m;
const spacesS = (props: CommonThemeProps) => getSpaces(props)?.s;
const spacesXxxs = (props: CommonThemeProps) => getSpaces(props)?.xxxs;

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
    font-weight: ${semiBold};
    margin-bottom: ${spacesS}px;
    padding: ${spacesXxxs}px 0;
  }

  td {
    padding: ${spacesM}px 0;
  }
`;
