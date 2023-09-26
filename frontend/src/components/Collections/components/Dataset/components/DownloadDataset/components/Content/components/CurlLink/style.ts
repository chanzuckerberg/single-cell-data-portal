import styled from "@emotion/styled";
import { fontBodyS, fontBodyXs } from "@czi-sds/components";
import {
  cornersM,
  grey100,
  grey300,
  grey500,
  spacesM,
  spacesS,
  spacesXxs,
} from "src/common/theme";

export const CodeBlock = styled.div`
  align-items: flex-start;
  background-color: ${grey100};
  border-radius: ${cornersM}px;
  box-shadow: inset 0 0 0 0.5px ${grey300};
  display: flex;
  gap: ${spacesS}px;
  margin: 0;
  padding: ${spacesS}px ${spacesM}px;

  code {
    ${fontBodyS}
    background-color: transparent;
    box-sizing: content-box;
    flex: 1;
    max-height: 40px;
    overflow: hidden;
    padding: ${spacesXxs}px 0 0;
  }

  code:before,
  code:after {
    content: unset; /* overrides styles from layout.css. */
  }
`;

export const Caption = styled.div`
  ${fontBodyXs}
  color: ${grey500};
`;
