import { Button } from "@czi-sds/components";
import styled from "@emotion/styled";

// We need this because this version of SDS does not respect the isAllCaps=false prop.
export const StyledButton = styled(Button)`
  text-transform: none;

  /* Apply SDS body-s-500 typography */
  font-size: 14px;
  font-weight: 500;
  letter-spacing: 0px;
  line-height: 24px;
  font-family:
    var(--font-inter),
    Inter,
    -apple-system,
    BlinkMacSystemFont,
    Segoe UI,
    Roboto,
    Helvetica Neue,
    Helvetica,
    Arial,
    sans-serif;
`;
