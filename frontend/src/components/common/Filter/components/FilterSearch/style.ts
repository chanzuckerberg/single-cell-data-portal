import { Classes, InputGroup } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { GRAY, LIGHT_GRAY, PT_TEXT_COLOR } from "src/components/common/theme";
import { getColors, IconButton } from "czifui";

export const ViewSearch = styled(InputGroup)`
  &.${Classes.INPUT_GROUP} {
    margin-bottom: 4px;

    .${Classes.ICON} {
      color: ${GRAY.A};
      margin: 0;
      padding: 4px;
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
      right: 8px;
    }
  }
`;

export const ClearIconButton = styled(IconButton)`
  ${(props) => {
    const colors = getColors(props);
    return `
      color: ${colors?.gray[400]};
      `;
  }};

  &:hover {
    color: black;
  }

  .MuiSvgIcon-root {
    height: 16px; // overrides icon size specification - xMarkCircle IconNameToSizes maximum available sdsSize is "s"
    width: 16px; // overrides icon size specification - xMarkCircle IconNameToSizes maximum available sdsSize is "s"
  }
`;
