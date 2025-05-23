import styled from "@emotion/styled";
import { fontCodeS } from "@czi-sds/components";
import {
  cornersM,
  grey100,
  grey300,
  spacesM,
  spacesS,
  spacesXxs,
} from "src/common/theme";
import { OLD_GRAY } from "src/components/common/theme";

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
    ${fontCodeS}
    background-color: transparent;
    box-sizing: content-box;
    flex: 1;
    line-height: 20px;
    overflow: hidden;
    padding: ${spacesXxs}px 0 0;
    white-space: pre-line;
  }

  code:before,
  code:after {
    content: unset; /* overrides styles from layout.css. */
  }
`;

export const DownloadCodeBlock = styled.div`
  margin-top: 14px;
  position: relative;
  width: 511px;

  code {
    border: 1px solid ${OLD_GRAY.VERY_LIGHT};
    border-radius: ${cornersM}px;
    box-sizing: border-box;
    display: block;
    font-size: 13px;
    letter-spacing: normal;
    line-height: 18px;
    max-height: 52px;
    overflow: hidden;
    padding: ${spacesS}px;
  }
`;
