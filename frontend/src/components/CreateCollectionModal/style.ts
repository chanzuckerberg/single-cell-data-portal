import styled from "@emotion/styled";
import { Button } from "@czi-sds/components";

export const StyledButton = styled(Button)`
  font-weight: 500;
  height: 32px; /* overrides @czi-sds/components height style declaration */
  justify-self: flex-end;
  letter-spacing: -0.006em;
  padding: 6px 12px;
`;
