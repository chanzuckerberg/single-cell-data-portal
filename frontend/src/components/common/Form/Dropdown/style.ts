import styled from "@emotion/styled";
import {
  CommonThemeProps,
  DropdownPopper as SDSDropdownPopper,
  getColors,
} from "@czi-sds/components";
import { css } from "@emotion/react";

interface DropdownFormProps extends CommonThemeProps {
  isSelected: boolean;
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

export const DropdownForm = styled("div")<DropdownFormProps>`
  .MuiButton-root {
    ${(props) => {
      const colors = getColors(props);
      return css`
        border-color: ${colors?.gray[400]};

        :hover {
          border-color: ${colors?.gray[400]};
        }
      `;
    }}

    height: 32px;
    width: 100%;
  }

  .MuiButton-text {
    > span {
      letter-spacing: -0.003em;

      ${(props) => {
        const textPrimary = props.theme?.palette?.text?.primary;
        return (
          props.isSelected &&
          css`
            color: ${textPrimary};
          `
        );
      }}
    }

    svg {
      height: 16px;
      width: 16px;

      path {
        ${(props) => {
          const colors = getColors(props);
          return (
            !props.isSelected &&
            css`
              fill: ${colors?.gray[500]};
            `
          );
        }}
      }
    }
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
