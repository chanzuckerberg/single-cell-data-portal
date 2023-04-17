import { Button, fontBodyS } from "czifui";
import styled from "@emotion/styled";

export const CollectionActions = styled.div`
  align-items: center;
  column-gap: 8px;
  display: flex;
`;

export const ActionButton = styled(Button)`
  ${fontBodyS};
  font-weight: 500;
  height: 32px;
  letter-spacing: -0.006em;
  min-width: 0;
  padding: 6px 12px;
`;
