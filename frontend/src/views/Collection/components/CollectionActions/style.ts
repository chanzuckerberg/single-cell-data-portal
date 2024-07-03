import styled from "@emotion/styled";
import { spacesS } from "src/common/theme";
import { Button } from "src/components/common/Button";

export const CollectionActions = styled.div`
  align-items: center;
  column-gap: ${spacesS}px;
  display: flex;
`;

export const ActionButton = styled(Button)`
  min-width: 0;
`;
