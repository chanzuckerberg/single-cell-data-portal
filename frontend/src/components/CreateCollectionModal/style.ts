import styled from "@emotion/styled";
import { Button } from "czifui";

export const StyledButton = styled(Button)`
  font-weight: 500;
  letter-spacing: -0.006em;
  height: 32px; /* overrides czifui height style declaration */
  justify-self: flex-end;
  padding: 6px 12px;
`;
