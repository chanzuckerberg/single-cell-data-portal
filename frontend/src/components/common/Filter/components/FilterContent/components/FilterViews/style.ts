import styled from "@emotion/styled";
import { CommonThemeProps, getSpaces } from "czifui";

const spacesXs = (props: CommonThemeProps) => getSpaces(props)?.xs;

export const ViewsMenu = styled.span`
  display: grid;
  flex: 1;
  grid-auto-columns: 1fr auto;
  grid-auto-flow: column;
  padding: ${spacesXs}px;
`;
