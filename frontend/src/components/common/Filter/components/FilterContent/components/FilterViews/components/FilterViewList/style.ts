import styled from "@emotion/styled";
import { List } from "src/components/common/Filter/components/FilterContent/components/common/style";
import { CommonThemeProps, getSpaces } from "@czi-sds/components";

const spacesXl = (props: CommonThemeProps) => getSpaces(props)?.xl;

export const ViewSublist = styled(List)`
  margin-left: ${spacesXl}px;
`;
