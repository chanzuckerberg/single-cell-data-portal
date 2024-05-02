import styled from "@emotion/styled";
import { fontBodyXs } from "@czi-sds/components";
import { Popper, PopperProps } from "@mui/material";

export const SectionItem = styled.li`
  ${fontBodyXs}
  font-weight: 400;
  cursor: pointer;
  padding-left: 20px !important;
  margin-bottom: 8px !important;
`;

interface StyledPopperProps extends PopperProps {
  fullWidth?: boolean;
}

export const StyledPopper = styled(Popper)<StyledPopperProps>`
  width: ${({ fullWidth }) => (fullWidth ? "100%" : "auto")};
  background-color: white;
  border: 1px solid #f8f8f8;
  border-radius: 4px;
  box-shadow:
    0 2px 4px 0 rgba(0, 0, 0, 0.15),
    0 2px 10px 0 rgba(0, 0, 0, 0.15);
  padding: 16px 4px 16px 16px;
  box-sizing: border-box;
  z-index: 1400;

  > .SdsAutocompleteMultiColumn-wrapper {
    height: 310px;

    > .SdsAutocompleteMultiColumn-column-root .base-Popper-root {
        transform: translate3d(0px, 20px, 0px) !important;
      }
    }
  }
`;
