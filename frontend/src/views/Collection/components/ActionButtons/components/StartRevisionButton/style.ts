import { Button, fontBodyS } from "czifui";
import styled from "@emotion/styled";

// TODO(cc) make generic button styles for action buttons.
export const StartRevisionButton = styled(Button)`
  ${fontBodyS};
  font-weight: 500;
  height: 32px;
  letter-spacing: -0.006em;
  min-width: 0;
  padding: 6px 12px;
`;
