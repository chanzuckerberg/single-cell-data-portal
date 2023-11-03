import styled from "@emotion/styled";
import { fontBodyS } from "@czi-sds/components";
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

/**
 * @deprecated by CodeBlock styles once "DOWNLOAD_UX" feature flag is removed.
 */
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
