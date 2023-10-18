import { Classes, InputGroup } from "@blueprintjs/core";
import { css } from "@emotion/react";
import styled from "@emotion/styled";
import { ButtonIcon, getColors } from "@czi-sds/components";
import { GRAY, LIGHT_GRAY, PT_TEXT_COLOR } from "src/components/common/theme";
import { spacesS, spacesXxs } from "src/common/theme";

export const ViewSearch = styled(InputGroup)`
  &.${Classes.INPUT_GROUP} {
    margin-bottom: ${spacesXxs}px;

    .${Classes.ICON} {
      color: ${GRAY.A};
      margin: 0;
      padding: ${spacesXxs}px;
      top: 50%;
      transform: translateY(-50%);
    }

    .${Classes.INPUT} {
      border: 1px solid ${LIGHT_GRAY.A} !important; /* required; overrides BP input border with important style declaration */
      border-radius: 3px;
      color: ${PT_TEXT_COLOR};
      letter-spacing: -0.1px;
      line-height: 18px;
      height: 32px;

      &::placeholder {
        opacity: 0.6;
      }
    }

    .${Classes.INPUT_ACTION} {
      right: ${spacesS}px;
      top: 50%;
      transform: translateY(-50%);
    }
  }
`;

export const ClearButtonIcon = styled(ButtonIcon)`
  ${(props) => {
    const colors = getColors(props);
    const palette = props.theme.palette;
    return css`
      color: ${colors?.gray[400]};

      &:hover {
        color: ${palette?.text?.primary};
      }
    `;
  }}
  .MuiSvgIcon-root {
    height: 16px; /* overrides IconNameToSizes "s" size specification */
    width: 16px; /* overrides IconNameToSizes "s" size specification */
  }
`;
