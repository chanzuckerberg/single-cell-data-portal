import styled from "@emotion/styled";
import { fontBodyXs } from "@czi-sds/components";

export const SectionItem = styled.li`
  ${fontBodyXs}
  font-weight: 400;
  cursor: pointer;
  padding-left: 20px !important;
  margin-bottom: 8px !important;
`;

interface WrapperProps {
  fullWidth?: boolean;
}
export const Wrapper = styled.div<WrapperProps>`
  .base-Popper-root {
    width: ${({ fullWidth }) => (fullWidth ? "100%" : "auto")};
    height: fit-content;
    box-shadow: 0 4px 4px 0 rgba(0, 0, 0, 0.25);
    border: none;
  }
  .base-Popper-root.MuiAutocomplete-popper.MuiAutocomplete-popperDisablePortal {
    box-shadow: none;
  }
`;
