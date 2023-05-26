import styled from "@emotion/styled";
import { ALIGNMENT } from "src/components/common/Grid/common/entities";
import { CommonThemeProps, getSpaces } from "@czi-sds/components";

const spacesXxs = (props: CommonThemeProps) => getSpaces(props)?.xxs;

interface Props extends CommonThemeProps {
  alignment: ALIGNMENT;
}

export const Header = styled("span")`
  align-items: center;
  align-self: center;
  display: flex;
  gap: ${spacesXxs}px; /* gap between header and sort icon */

  span {
    min-width: 0; /* facilitates breaking of word on columns; flex default for min width is "auto" */
  }
`;

export const HeaderCell = styled("th")<Props>`
  display: flex;
  gap: ${spacesXxs}px; /* gap between header and count */
  justify-content: ${(props) =>
    props.alignment === ALIGNMENT.LEFT ? "flex-start" : "flex-end"};
`;
