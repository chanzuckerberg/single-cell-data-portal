import styled from "@emotion/styled";
import { Popover } from "@mui/material";
import { shadowM, spacesXxs } from "src/common/theme";

export const Filter = styled.div`
  font-feature-settings: normal; /* required; overrides layout.css specification */
`;

export const FilterPopover = styled(Popover)`
  & .MuiPopover-paper {
    box-shadow: ${shadowM};
    font-feature-settings: normal; /* required; overrides layout.css specification */
  }
`;

export const InfoButtonWrapper = styled.span`
  cursor: pointer;
  margin-left: ${spacesXxs}px;
  /* width: 100% is needed to ensure its next/Image child has enough space to render the width and height specified */
  width: 100%;
`;
