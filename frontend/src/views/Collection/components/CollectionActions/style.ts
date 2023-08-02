import { Button, fontBodyS } from "@czi-sds/components";
import styled from "@emotion/styled";
import { spacesM, spacesS, spacesXs } from "src/common/theme";

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
