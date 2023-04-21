import styled from "@emotion/styled";
import { css } from "@emotion/react";
import { DropdownPopper as SDSDropdownPopper } from "czifui";

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

export const DropdownPopper = styled(SDSDropdownPopper)`
  border: none;

  .MuiAutocomplete-popper .MuiAutocomplete-paper {
    & .MuiAutocomplete-listbox {
      max-height: 224px; /* Displays max 7 options at 32px height each */
      padding-right: 4px; /* Allowance for scrollbar */

      .MuiAutocomplete-option {
        margin: 0;

        .MuiMenuItem-root {
          padding: 6px 8px;
          text-overflow: ellipsis;

          .MuiSvgIcon-root {
            align-self: center;
            margin: 0 10px 0 0;
          }

          & span {
            min-width: 0;

            div {
              letter-spacing: -0.003em;
              overflow: hidden;
              text-overflow: ellipsis;
              white-space: nowrap;
            }
          }
        }
      }

      ${scrollbar}
    }
  }
`;
