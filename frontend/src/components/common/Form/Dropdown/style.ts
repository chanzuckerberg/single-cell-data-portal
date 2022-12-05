import styled from "@emotion/styled";
import { DropdownPopper as SDSDropdownPopper } from "czifui";
import { css } from "@emotion/react";

const MENU_ITEM_HEIGHT_PX = 32;
export const MAX_DISPLAYABLE_MENU_ITEMS = 6;
const MAX_MENU_HEIGHT_PX = MAX_DISPLAYABLE_MENU_ITEMS * MENU_ITEM_HEIGHT_PX;

interface Props {
  isScrollable: boolean;
}

const scrollbar = css`
  &::-webkit-scrollbar {
    width: 8px;
  }

  &::-webkit-scrollbar-thumb {
    background-clip: content-box;
    background-color: rgba(0, 0, 0, 0.15);
    border-radius: 10px;
  }
`;

export const DropdownForm = styled.div`
  .MuiButton-root {
    height: 32px;
    margin: 0;
    width: 100%;
  }

  .MuiButton-label {
    min-width: 0; /* prevents button label in text overflow ellipsis from consuming full width of button container (flex default for min width is "auto") */

    svg {
      height: 16px;
      width: 16px;
    }
  }
`;

export const DropdownPopper = styled(SDSDropdownPopper, {
  shouldForwardProp: (prop) => prop !== "isScrollable",
})<Props>`
  border: none;

  .MuiAutocomplete-listbox {
    max-height: ${MAX_MENU_HEIGHT_PX}px;
    ${({ isScrollable }) => {
      return (
        isScrollable &&
        css`
          padding-right: 4px;
        `
      );
    }}
    ${scrollbar};
  }

  && .MuiAutocomplete-option {
    margin: 0;
  }

  .MuiListItem-root {
    letter-spacing: -0.003em;
    padding: 6px 8px;
    text-overflow: ellipsis;

    .MuiSvgIcon-root {
      margin: 0 10px 0 0;
    }

    & > span {
      min-width: 0;

      span {
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
      }
    }
`;
