import styled from "@emotion/styled";
import { CommonThemeProps, getColors, getSpaces } from "czifui";

const gray300 = (props: CommonThemeProps) => getColors(props)?.gray[300];
const spacesL = (props: CommonThemeProps) => getSpaces(props)?.l;

export const FilterDivider = styled.hr`
  background: none;
  box-shadow: inset 0px -0.5px 0px ${gray300};
  height: 0.5px;
  margin: ${spacesL}px 0;
`;
