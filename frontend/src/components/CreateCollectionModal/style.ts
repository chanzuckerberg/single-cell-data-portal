import styled from "@emotion/styled";
import { Button } from "czifui";

export const StyledButton = styled(Button)`
  font-weight: 500;
  height: 32px; /* overrides czifui height style declaration */
  justify-self: flex-end;
  letter-spacing: -0.006em;
  padding: 6px 12px;
`;
