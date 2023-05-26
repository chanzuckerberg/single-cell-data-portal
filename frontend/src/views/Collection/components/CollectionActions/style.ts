import { Button, CommonThemeProps, fontBodyS, getSpaces } from "czifui";
import styled from "@emotion/styled";

const spacesM = (props: CommonThemeProps) => getSpaces(props)?.m;
const spacesS = (props: CommonThemeProps) => getSpaces(props)?.s;
const spacesXs = (props: CommonThemeProps) => getSpaces(props)?.xs;

export const CollectionActions = styled.div`
  align-items: center;
  column-gap: ${spacesS}px;
  display: flex;
`;

export const ActionButton = styled(Button)`
  ${fontBodyS}
  font-weight: 500;
  height: 32px;
  letter-spacing: -0.006em;
  min-width: 0;
  padding: ${spacesXs}px ${spacesM}px;
`;
