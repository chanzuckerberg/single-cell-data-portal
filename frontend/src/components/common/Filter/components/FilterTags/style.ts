import { Classes, Colors } from "@blueprintjs/core";
import styled from "styled-components";

export const SelectedTags = styled.span`
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  min-width: 0; /* facilitates ellipsis on tag should it be required; flex default for min width is "auto" */

  .${Classes.TAG}.${Classes.MINIMAL} {
    background-color: #0073ff;
    border-radius: 4px;
    color: ${Colors.WHITE};
    font-size: 13px;
    font-weight: 500;
    letter-spacing: -0.1px;
    line-height: 18px;
    padding: 4px 8px;

    .${Classes.TAG_REMOVE} {
      margin: 0;
      opacity: 1;
      padding: 0;

      &:focus {
        outline: none;
      }

      .${Classes.ICON} {
        color: inherit;

        svg {
          height: 16px;
          width: 16px;
        }
      }
    }
  }

  .${Classes.TAG}.${Classes.LARGE} > * {
    margin-right: 8px;
  }
`;
