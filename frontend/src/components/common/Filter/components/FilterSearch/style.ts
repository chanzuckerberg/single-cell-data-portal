// TODO(cc) shared styles with FilterViewSearch - refactor or reuse FilterSearch.
import { Classes, InputGroup } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { GRAY, LIGHT_GRAY, PT_TEXT_COLOR } from "src/components/common/theme";

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
  }
`;
